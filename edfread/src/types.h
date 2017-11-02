#ifndef _E2M_TYPES
#define _E2M_TYPES

#include <vector>
#include <map>
#include <list>
#include <string>
using namespace std;

#ifndef MEX

#define mxClassID int
#define mxUINT32_CLASS 13
#define mxSINGLE_CLASS 7
#define mxINT32_CLASS 12

#endif /* MEX */

typedef map<string, string> strstrmap;

template<mxClassID C, typename T>
class StructField : public vector< T > {
	public:
		StructField( int size ) : vector< T >(size){ }
		StructField( ) : vector< T >( ){ }
		#ifdef MEX
		inline mxArray * toArray(T start=0) {
			int len = this->size();
			if (len==0)
				return mxCreateLogicalScalar( false );
			mxArray * result = mxCreateNumericMatrix(1, len, C, mxREAL); 
			T * out = (T*)mxGetData(result); 
			for (int i=0; i<len; i++){
				out[i] = (*this)[i]-start;
                        }
                        //cout << (*this)[0]-start;

			return result; 
		}
		#endif
};

class StringField : public vector< string > {
	public:
		#ifdef MEX
		inline mxArray * toArray() {
			const char ** msgstrs = (const char**)mxCalloc( this->size(), sizeof(char*) );
			for (int i=0; i<(int)this->size(); i++)
				msgstrs[i] = (*this)[i].c_str();
			mxArray * res = mxCreateCharMatrixFromStrings( this->size(), msgstrs );
			mxFree(msgstrs);
			return res;
		}
		#endif
};

typedef StructField<mxUINT32_CLASS,UINT32> UINT32Field;
typedef StructField<mxINT32_CLASS,INT32> INT32Field;
typedef StructField<mxSINGLE_CLASS,float> FloatField;

class Message {
	public:
		INT32Field time;
		StringField msg;
		#ifdef MEX
		inline mxArray * toArray(UINT32 start) {
			int dims[1] = { 1 };
			mxArray * res = mxCreateStructArray( 1, dims, MSG_COUNTF, msg_fields);
			mxSetFieldByNumber( res, 0, MSG_S, time.toArray(start) );
			mxSetFieldByNumber( res, 0, MSG_M, msg.toArray() );
			return res;
		}
		#endif
};
/*
 * Maps string to integer-string vector
 */
typedef map<string, Message> msgmap;

class EyeCalib {
	public:
		EyeCalib () : coeff(10) {
		avg =
		max =
		off_deg =
		off_y =
		res_x =
		res_y =
		off_x = 
#ifdef MEX
		mxGetNaN();
#else
		0.0;
#endif
		ctype = "ERR" ;
		}
		float avg, max, off_deg, off_x, off_y, res_x, res_y;
		string ctype;
		FloatField coeff;
#ifdef MEX
		inline mxArray * toArray() {
			int dims[1] = { 1 };
			mxArray * res = mxCreateStructArray( 1, dims, EC_COUNTF, ec_fields);
			mxSetFieldByNumber( res, 0, EC_AVG, mxCreateDoubleScalar( avg ) );
			mxSetFieldByNumber( res, 0, EC_MAX, mxCreateDoubleScalar( max ) );
			mxSetFieldByNumber( res, 0, EC_ODEG, mxCreateDoubleScalar( off_deg ) );
			mxSetFieldByNumber( res, 0, EC_OX, mxCreateDoubleScalar( off_x ) );
			mxSetFieldByNumber( res, 0, EC_OY, mxCreateDoubleScalar( off_y ) );
			mxSetFieldByNumber( res, 0, EC_RX, mxCreateDoubleScalar( res_x ) );
			mxSetFieldByNumber( res, 0, EC_RY, mxCreateDoubleScalar( res_y ) );
			const char * stype = ctype.c_str();
			mxSetFieldByNumber( res, 0, EC_TYPE, mxCreateCharMatrixFromStrings(1, &stype ) );
			mxSetFieldByNumber( res, 0, EC_COEFF, coeff.toArray() );
			return res;
		}		
#endif
};

class Calibration {
	public:
		vector<EyeCalib *> eyes;
		int btrial;
		Calibration( int trial ) : eyes(2) {
			eyes[0] = NULL;
			eyes[1] = NULL;
			btrial = trial;
		}
		~Calibration(){
			if (eyes[0]!=NULL) delete eyes[0];
			if (eyes[1]!=NULL) delete eyes[1];
		}
		bool add_eye( int id ) {
			if (eyes[id] == NULL) {
				eyes[id] = new EyeCalib();
				return true;
			} else return false;
		}
		EyeCalib * get_eye( int id ) {
			if (eyes[id] == NULL) 
				eyes[id] = new EyeCalib();
			return eyes[id];
		}
};

typedef list<Calibration *> calib_list; 

class Button {
	public:
		INT32Field time;
		UINT32Field code;
		#ifdef MEX
		inline mxArray * toArray(UINT32 start) {
			int dims[1] = { 1 };
			mxArray * res = mxCreateStructArray( 1, dims, BTN_COUNTF, btn_fields);
			mxSetFieldByNumber( res, 0, BTN_T, time.toArray(start) );
			mxSetFieldByNumber( res, 0, BTN_K, code.toArray() );
			return res;
		}
		#endif
};


class Fixation {
	public:
		INT32Field start;
		INT32Field end;
		FloatField x;
		FloatField y;
		FloatField pupil;
		#ifdef MEX
		inline mxArray * toArray(UINT32 tstart) {
			int dims[1] = { 1 };
			mxArray * res = mxCreateStructArray( 1, dims, FIX_COUNTF, fix_fields);
			mxSetFieldByNumber( res, 0, FIX_S, start.toArray(tstart) );
			mxSetFieldByNumber( res, 0, FIX_E, end.toArray(tstart) );
			mxSetFieldByNumber( res, 0, FIX_X, x.toArray() );
			mxSetFieldByNumber( res, 0, FIX_Y, y.toArray() );
			mxSetFieldByNumber( res, 0, FIX_P, pupil.toArray() );
			return res;
		}
		#endif
};

class Saccade {
	public:
		INT32Field start;
		INT32Field end;
		FloatField sx;
		FloatField ex;
		FloatField sy;
		FloatField ey;
		FloatField speed;
		#ifdef MEX
		inline mxArray * toArray(UINT32 tstart) {
			int dims[1] = { 1 };
			mxArray * res = mxCreateStructArray( 1, dims, SAC_COUNTF, sac_fields);
			mxSetFieldByNumber( res, 0, SAC_ST, start.toArray(tstart) );
			mxSetFieldByNumber( res, 0, SAC_ET, end.toArray(tstart) );
			mxSetFieldByNumber( res, 0, SAC_SX, sx.toArray() );
			mxSetFieldByNumber( res, 0, SAC_SY, sy.toArray() );
			mxSetFieldByNumber( res, 0, SAC_EX, ex.toArray() );
			mxSetFieldByNumber( res, 0, SAC_EY, ey.toArray() );
			mxSetFieldByNumber( res, 0, SAC_SP, speed.toArray() );
			return res;
		}
		#endif
};

class Blink {
	public:
		INT32Field start;
		INT32Field end;
		#ifdef MEX
		inline mxArray * toArray(UINT32 tstart) {
			int dims[1] = { 1 };
			mxArray * res = mxCreateStructArray( 1, dims, BLINK_COUNTF, blink_fields);
			mxSetFieldByNumber( res, 0, BLINK_START, start.toArray(tstart) );
			mxSetFieldByNumber( res, 0, BLINK_END, end.toArray(tstart) );
			return res;
		}
		#endif
};

class Sample {
	public:
		INT32Field time;
		FloatField x;
		FloatField y;
		FloatField pupil;
		#ifdef MEX
		inline mxArray * toArray(UINT32 start) {
			int dims[1] = { 1 };
			mxArray * res = mxCreateStructArray( 1, dims, SAMPLE_COUNTF, sample_fields);
			mxSetFieldByNumber( res, 0, SAMPLE_TIME, time.toArray(start) );
			mxSetFieldByNumber( res, 0, SAMPLE_X, x.toArray() );
			mxSetFieldByNumber( res, 0, SAMPLE_Y, y.toArray() );
			mxSetFieldByNumber( res, 0, SAMPLE_PUPIL, pupil.toArray() );
			return res;
		}
		#endif
};

#endif /* _E2M_TYPES */
