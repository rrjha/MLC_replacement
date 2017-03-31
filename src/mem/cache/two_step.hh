#ifndef __TWO_STEP_H__
#define __TWO_STEP_H__

#include "base/statistics.hh"

using namespace std;

/********************** Transition Energy ***********************
From/To |   R00     R01     R10     R11
--------------------------------------------
R00     |   ZT      ST      TT      HT
R01     |   ST      ZT      TT      HT
R10     |   HT      TT      ZT      ST
R11     |   HT      TT      ST      ZT
*********************************************************************/
enum {ZT, ST, HT, TT, MAX_TRANSITION};

typedef struct decision_table_entry {
    uint8_t code;
    uint32_t transitions[MAX_TRANSITION];
} dtab_entry;

class two_step
{
    public:
        two_step(Stats::Vector& trans);
        virtual ~two_step();
        void write_ts_encoded(uint8_t *, const uint8_t *, uint32_t, bool ignore_energy = false);

    private:

        void find_energy_score(int8_t, int8_t, uint32_t *);
        void gen_table();
        void decode (uint32_t, uint16_t *);
        void encode(uint16_t, uint32_t *, bool ignore_energy = false);
        void read_ts_decoded(const uint8_t *, uint8_t *, uint32_t);

        const float energy[MAX_TRANSITION] = { /* ZT */ 0, /* ST */ 1.084, /* HT */ 2.653, /* TT */ 3.737 };
        const int8_t codetab[4][3] =	{ {   0, -1, -1   }, {   1,  2,  4   }, {   3,  5,  6   }, {   7, -1, -1   }   };
        static dtab_entry m_decisiontab[64][16];
        Stats::Vector &m_transitions; //assume the allocation is done for all transitions by caller
};

#endif // __TWO_STEP_H__
