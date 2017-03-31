#include "two_step.hh"
#include "base.hh"
dtab_entry two_step::m_decisiontab[64][16];
//Stats::Scalar* two_step::m_transitions;

two_step::two_step(Stats::Vector& trans):m_transitions(trans)
{


    gen_table();
}

two_step::~two_step()
{
    //dtor
    //m_transitions = NULL;
}


void two_step::find_energy_score(int8_t from_code, int8_t to_code, uint32_t *transitions) {
    uint8_t i, from_bits, to_bits;
    for(i=0; i < 6; i+=2) {
        from_bits = (from_code >> i) & 3;
        to_bits = ((to_code >> i) & 3);
        switch(from_bits) {
            case 0:
                (to_bits == 1) ? transitions[ST]++ :
                (to_bits == 2) ? transitions[TT]++ :
                (to_bits == 3) ? transitions[HT]++ : transitions[ZT]++;
                break;
            case 1:
                (to_bits == 2) ? transitions[TT]++ :
                (to_bits == 3) ? transitions[HT]++ :
                (to_bits == 0) ? transitions[ST]++ : transitions[ZT]++;
                break;
            case 2:
                (to_bits == 3) ? transitions[ST]++ :
                (to_bits == 0) ? transitions[HT]++ :
                (to_bits == 1) ? transitions[TT]++ : transitions[ZT]++;
                break;
            case 3:
                (to_bits == 0) ? transitions[HT]++ :
                (to_bits == 1) ? transitions[TT]++ :
                (to_bits == 2) ? transitions[ST]++ : transitions[ZT]++;
                break;
        }
    }
}

void two_step::gen_table() {
    uint8_t i, j, lo, hi, m, n, p, code;
    uint32_t transitions[MAX_TRANSITION];
    float energy_score, least_energy_score;
    dtab_entry chosen;

    for(i=0; i<64; i++) {
        for(j=0; j<16; j++) {
            /* For each incoming 4bit dissect in two and generate 3 bit possible patterns
             * for a total of max 9 6-bit patterns for this 4 bit */
             lo = j & 3;
             hi = ((j>>2) & 3);
             code = 0;
             memset(&chosen, 0, sizeof(dtab_entry));
             least_energy_score = 3*energy[TT]; //max possible score to start with

             for(m=0; (m < 3)&&(codetab[hi][m] != -1); m++) {
                for(n=0; (n < 3)&&(codetab[lo][n] != -1); n++) {
                    code = ((codetab[hi][m] << 3) | (codetab[lo][n])) & 0x3F;
                    memset(transitions, 0, sizeof(uint32_t)*MAX_TRANSITION);
                    find_energy_score(i, code, transitions);
                    energy_score = transitions[ZT] * energy[ZT] + transitions[ST] * energy[ST] +
                                    transitions[HT] * energy[HT] + transitions[TT] * energy[TT];
                    if(energy_score < least_energy_score) {
                        least_energy_score = energy_score;
                        chosen.code = code;
                        for(p=0; p < MAX_TRANSITION; p++)
                            chosen.transitions[p] = transitions[p];
                    }
                }
            }
            m_decisiontab[i][j] = chosen;
        }
    }
}

void two_step::decode (uint32_t fromdata, uint16_t *p_todata) {
    int i=0;
    uint16_t result=0;
    uint8_t curr_three_bits = 0;
    uint8_t set_bits = 0;

    for(i=0; i<8; i++) { // Loop in 8*3 bits
        curr_three_bits = (fromdata >> ((7 - i) * 3)) & 7; //right shift and AND with 111 to get the three bits for this iter
        set_bits = 0;
        while(curr_three_bits) { // find number of set bits for decoding this three bits to two actual bits
            set_bits++;
            curr_three_bits &= (curr_three_bits-1);
        }
        result = (result << 2) | set_bits; //append 2 bits to result
    }
    *p_todata = result;
}


void two_step::encode(uint16_t todata, uint32_t *p_fromdata, bool ignore_energy){
    uint32_t i=0, j=0;
    uint32_t result = (*p_fromdata) & 0xFF000000;; //retain MSB
    dtab_entry temp;

    for(i=0; i<4; i++) { //loop for 4 4-bit nibbles converting them using decision table
        temp = m_decisiontab[((*p_fromdata >> 6*i) & 0x3F)][((todata >> 4*i) & 0xF)];
        result |= (temp.code << 6*i);
        if(!ignore_energy) {
            for(j = 0; j < MAX_TRANSITION; j++)
                m_transitions[j] += temp.transitions[j];
        }
    }
    *p_fromdata = result;
}

void two_step::write_ts_encoded(uint8_t *fromblk, const uint8_t *toblk, uint32_t blksize, bool ignore_energy) {
    int32_t i=0, j=0;
    uint32_t residual = 0;

    for(i=0, j=0; i < (blksize-2); i+=2, j+=3)
        encode(*((uint16_t*)(toblk+i)), (uint32_t*)(fromblk+j), ignore_energy);
    residual |= ((*(fromblk+j+2) << 16) | (*(fromblk+j+1) << 8) | *(fromblk+j)); //last 3 bytes remain
    encode(*((uint16_t*)(toblk+i)), (uint32_t*)&residual, ignore_energy);
    std::memcpy((fromblk+j), &residual, 3);
}

void two_step::read_ts_decoded(const uint8_t *fromblk, uint8_t *toblk, uint32_t blksize) {
    uint32_t i=0, j=0;

    for(i=0, j=0; i < (blksize + (blksize >>1)); i+=3, j+=2) //Size of fromdata is 3/2 times of target so adjust size and then keep 3 space for pointer in loop
        decode(*((uint32_t*)(fromblk+i)) & 0x00FFFFFF, (uint16_t*)(toblk+j));
}


