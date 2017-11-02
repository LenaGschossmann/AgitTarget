#include "types.h"
#include <string.h>

/*
 * Represents all data for one Trial
 */
class Trial {
public:
    /* start time */
    UINT32 start;
    /* button events */
    Button ebuttons;
    /* map for storing filter results */
    msgmap msgs;
    /* vector of two Eyes */
    Eyes * eyes;
    /*
     * Initialize, filter messages for strings given in "fliter"
     */
    Trial(char * filter[], int flen) : ebuttons(), msgs()
    {
        eyes = new Eyes(2);
				start = 0;
        for (int i=0; i<flen; i++)
        {
            Message m;
            msgs.insert(msgmap::value_type(filter[i], m));
        }
    }
    ~Trial() {
        delete eyes;
    }
    /*
     * forward a new sample to the corresponding Eye
     */
    inline void sample(FSAMPLE * sam){
        if (sam->flags & SAMPLE_RIGHT)
            (*eyes)[1].sample(sam->time, sam->gx[1], sam->gy[1], sam->pa[1]);
        if (sam->flags & SAMPLE_LEFT)
            (*eyes)[0].sample(sam->time, sam->gx[0], sam->gy[0], sam->pa[0]);
    }
    /*
     * adds a button press
     */
    inline void button(UINT32 time, UINT32 button){
        ebuttons.time.push_back(time);
        ebuttons.code.push_back(button);
    }
    /*
     * adds a filtered message
     */
    inline void message(char * filter, UINT32 time, const char * msg){
        msgs[filter].time.push_back((unsigned int)time);
        msgs[filter].msg.push_back(string(msg));
    }
    /*
     * add a metadata message
     */
    inline void meta(string key, UINT32 time, string msg){
	if (msgs.find(key)==msgs.end())
	{
            Message m;
            msgs.insert(msgmap::value_type(key, m));
	}
        msgs[key].time.push_back(time);
	msgs[key].msg.push_back(msg);
    }

};

typedef list <Trial *> trial_list;
