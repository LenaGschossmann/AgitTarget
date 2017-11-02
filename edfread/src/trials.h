#include "edf.h"
#include <iostream>
#include <sstream>
#include "types.h"
#include "eye.h"
#include "trial.h"
#include "string.h"

#ifdef MEX
#define DEBUG(x)	{};
#define DEBUG2(x)	{};
#define DEBUGNSW(x)     {};
#else
#define DEBUG(x)	cerr << __LINE__ << " " << x << endl;
#define DEBUG2(x)	cerr << __LINE__ << ":2 " << x << endl;
#define DEBUGNSW(x)  cerr << __LINE__ << ":2 " << x << endl;

#endif

using namespace std;

/*
 * Main access class for EDF,
 * builds and maintains list of Trials
 * for parsing and mxArray retrieval
 */
class Trials
{
private:
    char **msgfilter;
    int filtermsgcount;
    Trial *current;
    trial_list *trials;
    calib_list calib;
    string header;
    int current_calibe;
    Calibration *current_calib;
    int start_count;
    strstrmap metadata;

public:
    Trials(char *filter[], int flen)
    {
        filtermsgcount = flen;
        msgfilter = filter;
        current = NULL;
        current_calib = NULL;
        trials = new trial_list;
        /* sengmann
        newTrial();
        cerr << "WARNING: special sonja version for old c experiments" << endl;
        */
        start_count = 0;
    }
    ~Trials()
    {
        current = NULL;
        while (trials->size() > 0)
        {
            Trial *t = trials->back();

            trials->pop_back();
            delete t;
        }
        while (calib.size() > 0)
        {
            Calibration *t = calib.back();

            calib.pop_back();
            delete t;
        }
        delete trials;
    }
#ifdef MEX
    mxArray *data_struct()
    {
        int dims[2] = {1, trials->size()};
        mxArray *result = mxCreateStructArray(2, dims, TRIAL_COUNTF, trial_fields);
        for (int i = 0; i < filtermsgcount; i++)
            mxAddField(result, (const char *) msgfilter[i]);
        int i = 0;
        for (trial_list::const_iterator it = trials->begin();
                it != trials->end();
                it++, i++)
        {
            Trial *t = (*it);
            mxSetFieldByNumber(result, i, TRIAL_LEFT, (*(t->eyes))[0].asstruct(t->start));
            mxSetFieldByNumber(result, i, TRIAL_RIGHT, (*(t->eyes))[1].asstruct(t->start));
            mxSetFieldByNumber(result, i, TRIAL_BUTTON, t->ebuttons.toArray(t->start));
            for (msgmap::iterator q = t->msgs.begin(); q != t->msgs.end(); q++)
            {
                mxAddField(result, (*q).first.c_str());
                mxSetField(result, i, (*q).first.c_str(), (*q).second.toArray(t->start));
            }
        }
        return result;
    }
    mxArray *info_struct()
    {
        int dims[1] = {1};
        mxArray *result = mxCreateStructArray(1, dims, INFO_COUNTF, info_fields);
        mxSetFieldByNumber(result, 0, INFO_HEAD, mxCreateString(header.c_str()));

        dims[0] = calib.size();
        mxArray *cal = mxCreateStructArray(1, dims, CAL_COUNTF, cal_fields);
        int i = 0;
        for (calib_list::const_iterator it = calib.begin();
                it != calib.end();
                it++, i++)
        {
            mxSetFieldByNumber(cal, i, CAL_TRIAL, mxCreateDoubleScalar((*it)->btrial));
            if ((*it)->eyes[0])
                mxSetFieldByNumber(cal, i, CAL_LEFT, (*it)->eyes[0]->toArray());
            else
                mxSetFieldByNumber(cal, i, CAL_LEFT, mxCreateLogicalScalar(false));
            if ((*it)->eyes[1])
                mxSetFieldByNumber(cal, i, CAL_RIGHT, (*it)->eyes[1]->toArray());
            else
                mxSetFieldByNumber(cal, i, CAL_RIGHT, mxCreateLogicalScalar(false));
        }

        mxSetFieldByNumber(result, 0, INFO_CAL, cal);

        for (strstrmap::iterator q = metadata.begin();
                q != metadata.end();
                q++)
        {
            mxAddField(result, (*q).first.c_str());
            mxSetField(result, 0, (*q).first.c_str(), mxCreateString((*q).second.c_str()));
        }
        return result;
    }

#endif

    void newCalibration()
    {
        current_calib = new Calibration(start_count);
        calib.push_back(current_calib);
    }

    void newTrial()
    {
        DEBUG("new trial " << start_count)
        current = new Trial(msgfilter, filtermsgcount);
        trials->push_back(current);
    }

    bool parse_meta(char *m, string & k, string & v)
    {
        string msg(m);
        msg.erase(0, 7);
        int f = msg.find(' ');
        if (f < 0)
            return false;
        k = msg.substr(0, f);
        v = msg.substr(f+1, msg.length()-f);
        DEBUG("METAEX " << k << "-" << v);
        return true;
    }

    void check_calib(FEVENT fe)
    {
        DEBUGNSW("IN CALIB");
        DEBUGNSW(&fe.message->c)
        /*
        * !CAL >>>>>>> CALIBRATION (HV9,P-CR) FOR LEFT: <<<<<<<<<
        */
        //if (!fe.message) return;
        if (strstr(&fe.message->c, ">>>>>>> CALIBRATION"))
        {
            DEBUGNSW("!strncmp(&fe.message->c, >>>>>>> CALIBRATION,20)")
            if (current_calib == NULL)
                newCalibration();
            istringstream s(&fe.message->c);
            string tmp;
            for (int i = 0; i < 6; i++)
                s >> tmp;
            current_calibe = tmp == "LEFT:" ? 0 : 1;
            if (!current_calib->add_eye(current_calibe))
            {
                DEBUG("New Calib")
                newCalibration();
                current_calib->add_eye(current_calibe);
            }
        }
        else
            DEBUGNSW("strncmp(&fe.message->c, >>>>>>> CALIBRATION,20)")
            /*
            * !CAL VALIDATION HV9 LR LEFT GOOD ERROR 0.47 avg.
            * 0.87 max OFFSET 0.27 deg. 4.2,-7.6 pix.
            */
            /* !CAL VALIDATION LR ABORTED */
            if (!strncmp(&fe.message->c, "!CAL VALIDATION", 15))
            {
                DEBUGNSW("CAL VALIDATION")
                istringstream s(&fe.message->c);
                string tmp;
                string ctype;

                s >> tmp;
                s >> tmp;
                s >> ctype;	/* type */
                s >> tmp;
                if (tmp != "ABORTED")
                {
                    s >> tmp;	/* eye */
                    int eye = tmp == "LEFT" ? 0 : 1;
                    if (current_calib != 0){
                        EyeCalib *ec = current_calib->get_eye(eye);
                        ec->ctype = ctype;
                        s >> tmp;
                        s >> tmp;
                        s >> ec->avg;	/* avg err */
                        s >> tmp;
                        s >> ec->max;	/* max arr */
                        s >> tmp;
                        s >> tmp;
                        s >> ec->off_deg;	/* offset */
                        s >> tmp;
                        char foo[10];

                        s.getline(foo, 10, ',');
                        ec->off_x = atof(foo);	/* pix x */
                        s >> ec->off_y;	/* pix y */
                    }else{
                        DEBUG("VALIDATION BUT NO CALIBRATION -> Is this correct?")
                    }
                }
                DEBUG("VAL DONE")
            }
            else
            {
                /*
                * !CAL Cal
                * coeff:(X=a+bx+cy+dxx+eyy,Y=f+gx+goaly+ixx+jyy)
                * 728.6 75.363 0.53358 0.093817 2.6942 1052.5
                * -1.1784 322.43 0.012614 3.095
                */
                DEBUGNSW("NOT VALIDATION")
                if (!strncmp(&fe.message->c, "!CAL Cal coeff", 14))
                {
                    DEBUGNSW("CAL COEFFs")
                    DEBUGNSW("PRINT MSG")
                    DEBUGNSW(&fe.message->c)
                    istringstream s(&fe.message->c);
                    string tmp;

                    s >> tmp;
                    s >> tmp;
                    s >> tmp;
                    for (int i = 0; i < 10; i++)
                        s >> current_calib->get_eye(current_calibe)->coeff[i];
                }
            }
        /* !CAL Resolution (upd) at screen center: X=3.5, Y=0.8 */
        if (!strncmp(&fe.message->c, "!CAL Resolution", 15))
        {
            EyeCalib *ec = current_calib->get_eye(current_calibe);
            istringstream s(&fe.message->c);
            string tmp;

            s >> tmp;
            s >> tmp;
            s >> tmp;
            s >> tmp;
            s >> tmp;
            s >> tmp;
            char foo[10];

            s.getline(foo, 10, '=');
            s.getline(foo, 10, ',');
            ec->res_x = atof(foo);	/* resolution x */
            s.getline(foo, 10, '=');
            s >> ec->res_y;	/* resolution y */
        }
    }

    void handle_message(FEVENT fe)
    {
        if (fe.message)
            if (fe.message->len > 0)
            {
                char *msg = &fe.message->c;
                DEBUG("Message" << msg);
                /* sengmann special*/
                if (!strncmp(msg, "TRIALID", 7))
                    /*
                    if (strstr(msg, "IMGLOAD"))*/
                {
                    DEBUG(msg);
                    newTrial();
                }
                else if (!strncmp(msg, "METAEX", 6))
                {
                    string k, v;
                    if (parse_meta(msg, k, v))
                        metadata[k] = v;
                }
                /* sengmann else */ if (current)
                {
                    if (!strncmp(msg, "SYNCTIME", 8))
                    {
                        if (current->start==0)
                            current->start = fe.sttime;
                    }
                    else if (!strncmp(msg, "DRIFTCORRECT", 12))
                    {
                        istringstream s(msg);
                        string tmp;

                        s >> tmp;
                        s >> tmp;
                        s >> tmp;
                        Eye *eye = &(*(current->eyes))[tmp == "LEFT" ? 0 : 1];

                        s >> tmp;
                        s >> tmp;
                        s >> tmp;
                        s >> eye->drift[2];
                        s >> tmp;
                        char foo[10];

                        s.getline(foo, 11, ',');
                        eye->drift[0] = atof(foo);
                        s >> eye->drift[1];
                        DEBUG("Drift " << eye->drift[0] << " " << eye->drift[1] << " " << eye->drift[2]);
                    }
                    else if (!strncmp(msg, "METATR", 6))
                    {
                        string k, v;
                        if (parse_meta(msg, k, v))
                            current->meta(k, fe.sttime, v);
                    }
                    else
                    {
                        for (int i = 0; i < filtermsgcount; i++)
                            if (strstr(msg, msgfilter[i]))
                            {
                                current->message(msgfilter[i], fe.sttime, msg);
                                break;
                            }
                    }
                }
            }
    }

    /*
    * Opens and parses an EDF into the instances structure
    */
    bool readfile(char *file)
    {
        EDFFILE *ed = NULL;
        bool active_trial = false;
        int err = 0;

        ed = edf_open_file(file, 0, 1, 1, &err);
        if (ed == NULL)
            return false;
        char buf[1024] = "";
        edf_get_preamble_text(ed, buf, 1024);
        header = buf;

        ALLF_DATA *fd = NULL;

        /* Main parsing loop */
        while (1)
        {
            int type = edf_get_next_data(ed);
            switch (type)
#define CUREYE ((*(current->eyes))[fd->fe.eye])
#define EVENT (fd->fe)
            {
            case STARTBLINK:
            case STARTSACC:
                break;
            case STARTFIX:
            {
                if (!current)
                    break;
                DEBUG2("StartFIX");
                fd = edf_get_float_data(ed);
                CUREYE.init_fix();
                break;
            }
            case STARTPARSE:
                break;

            case ENDBLINK:
            {
                if (!current)
                    break;
                DEBUG2("EndBLINK");
                fd = edf_get_float_data(ed);
                CUREYE.add_blink(EVENT.sttime, EVENT.entime);
                break;
            }
            case ENDSACC:
            {
                if (!current)
                    break;
                DEBUG2("EndSAC");
                fd = edf_get_float_data(ed);
                CUREYE.add_sac(EVENT.sttime,  EVENT.gstx, EVENT.gsty, EVENT.entime, EVENT.genx, EVENT.geny, EVENT.avel);
                break;
            }
            case ENDFIX:
            {
                if (!current)
                    break;
                DEBUG2("EndFix");
                fd = edf_get_float_data(ed);
                CUREYE.stop_fix(EVENT.sttime, EVENT.entime);
                break;
            }
            case ENDPARSE:
            case FIXUPDATE:
            case BREAKPARSE:
                break;
            case BUTTONEVENT:
            {
                DEBUG2("Button " << EVENT.sttime << " code " << EVENT.buttons);
                if (!current)
                    break;
                fd = edf_get_float_data(ed);
                int btn = EVENT.buttons >> 8;

                if (!btn)
                    break;
                int down = EVENT.buttons & 255;

                for (int i = 1; i < 9; i++)
                {
                    if (btn & 1)
                        if (down & 1)
                            current->button(EVENT.sttime, i);
                        else
                            current->button(EVENT.sttime, 10 + i);
                    btn >>= 1;
                    down >>= 1;
                }
                break;
            }
            case INPUTEVENT:
                break;
            case MESSAGEEVENT:
            {
                fd = edf_get_float_data(ed);
                if(!EVENT.message->c){
                    DEBUGNSW("NULLMSG");
                    return false;
                }
                handle_message(EVENT);
                DEBUGNSW("\n---------------------------------------")
                DEBUGNSW("Print Message: ")
                if (!active_trial)
                    check_calib(EVENT);
                break;
            }
            case STARTSAMPLES:
            case STARTEVENTS:
            case ENDSAMPLES:
                break;
            case ENDEVENTS:
                break;
            case SAMPLE_TYPE:
            {
                if (!current)
                    break;
                fd = edf_get_float_data(ed);
                current->sample(&(fd->fs));
                break;
            }
            case RECORDING_INFO:
                fd = edf_get_float_data(ed);
                active_trial = fd->rec.state == 1;
                DEBUG("RecordingInfo " << active_trial << " " << int(fd->rec.pos_type));
                if (active_trial)
                {
                    start_count++;
                }
                break;
            case NO_PENDING_ITEMS:
                DEBUG("NoPending");
                edf_close_file(ed);
                return true;
            default:
                /* a type of message we do not understand */
                edf_close_file(ed);
                return false;
            }
        }
    }

};
