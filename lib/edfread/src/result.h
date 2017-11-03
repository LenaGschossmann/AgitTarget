#define TRIAL_COUNTF	3
#define TRIAL_LEFT 	0
#define TRIAL_RIGHT	1
#define TRIAL_BUTTON	2
const char * trial_fields[TRIAL_COUNTF] = { "left", "right", "button" };

#define EYE_COUNTF	5
#define	EYE_FIXATION	0
#define	EYE_SACCADE	1
#define EYE_BLINK	2
#define EYE_SAMPLES	3
#define EYE_DRIFT	4
const char * eye_fields[EYE_COUNTF] = { "fixation", "saccade", "blink", "samples", "drift" };

#define FIX_COUNTF 5
#define FIX_S 0
#define FIX_E 1
#define FIX_X 2
#define FIX_Y 3
#define FIX_P 4
const char * fix_fields[FIX_COUNTF] = { "start", "end", "x", "y", "pupil" };

#define SAC_COUNTF 7
#define SAC_ST 0
#define SAC_SX 1
#define SAC_SY 2
#define SAC_ET 3
#define SAC_EX 4
#define SAC_EY 5
#define SAC_SP 6
const char * sac_fields[SAC_COUNTF] = { "start", "sx", "sy", "end", "ex", "ey", "speed" };

#define BLINK_COUNTF 2
#define BLINK_START 0
#define BLINK_END 1
const char * blink_fields[BLINK_COUNTF] = { "start", "end" };

#define SAMPLE_COUNTF 4
#define SAMPLE_TIME 0
#define SAMPLE_X 1
#define SAMPLE_Y 2
#define SAMPLE_PUPIL 3
const char * sample_fields[SAMPLE_COUNTF] = { "time", "x", "y", "pupil" };

#define MSG_COUNTF 2
#define MSG_S 0
#define MSG_M 1
const char * msg_fields[MSG_COUNTF] = { "time", "msg" };

#define BTN_COUNTF 2
#define BTN_T 0
#define BTN_K 1
const char * btn_fields[BTN_COUNTF] = { "time", "code" };

#define INFO_COUNTF 2
#define INFO_HEAD 0
#define INFO_CAL 1
const char * info_fields[INFO_COUNTF] = { "header", "calib" };

#define CAL_COUNTF 3
#define CAL_TRIAL 0
#define CAL_LEFT 1
#define CAL_RIGHT 2
const char * cal_fields[CAL_COUNTF] = { "trial", "left", "right" };

#define EC_COUNTF 9
#define EC_AVG 0
#define EC_MAX 1
#define EC_ODEG 2
#define EC_OX 3
#define EC_OY 4
#define EC_RX 5
#define EC_RY 6
#define EC_TYPE 7
#define EC_COEFF 8
const char * ec_fields[EC_COUNTF] = { "err_avg", "err_max", "off_deg", "off_x", "off_y", "res_x", "res_y", "type", "coeff" };
