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
        two_step(Stats::Scalar *);
        virtual ~two_step();
        static void write_ts_encoded(uint8_t *, const uint8_t *, uint32_t);

    private:

        void find_energy_score(int8_t, int8_t, uint32_t *);
        void gen_table();
        void decode (uint32_t, uint16_t *);
        static void encode(uint16_t, uint32_t *);
        void read_ts_decoded(const uint8_t *, uint8_t *, uint32_t);

        const float energy[MAX_TRANSITION] = { /* ZT */ 0, /* ST */ 1.92, /* HT */ 3.192, /* TT */ 5.112 };
        const int8_t codetab[4][3] =	{ {   0, -1, -1   }, {   1,  2,  4   }, {   3,  5,  6   }, {   7, -1, -1   }   };
        static dtab_entry m_decisiontab[64][16];
        static Stats::Scalar *m_transitions; //assume the allocation is done for all transitions by caller
};

#endif // __TWO_STEP_H__
