#include "types.h"

/*
 * Represents Eye-specific data for one Trial
 */
class Eye {
	public:
		Saccade sac;
		Fixation fix;
		Blink blink;
		Sample samples;
		FloatField drift;
		int start_fix;

		Eye() : drift(3) {
		start_fix = -1; 
		drift[0] = drift[1] = drift[2] = NAN;};

		/* add a saccade */
		inline void add_sac(UINT32 start, float sx, float sy, 
							UINT32 stop, float ex, float ey, float vel) { 
			sac.start.push_back( start ); 
			sac.sx.push_back( sx ); 
			sac.sy.push_back( sy ); 
			sac.end.push_back( stop ); 
			sac.ex.push_back( ex ); 
			sac.ey.push_back( ey ); 
			sac.speed.push_back( vel ); 
		} 

		/* add a blink */
		inline void add_blink(UINT32 start, UINT32 stop) {
			blink.start.push_back(start);
			blink.end.push_back(stop);
		}

		/* remember when the current fixation started */
		inline void init_fix() {
			start_fix = samples.x.size();
		}

		/*
		 * signal end of current fixation
		 * will calculate the fixation coordinate by
		 * triangularly weighting all samples during fixation
		 */
		inline void stop_fix(UINT32 start, UINT32 stop) {
		  if (start_fix==-1) return;
			const int ssize = samples.x.size();
			int len = ssize - start_fix;
			int mid = (int)trunc(len/2);
			double x = 0;
			double y = 0;
			double p = 0;
			for (int i=0; i<mid; i++) { 
				double w = (i+1.0)/mid;
				x += w*samples.x[start_fix+i];
				y += w*samples.y[start_fix+i];
				p += w*samples.pupil[start_fix+i];
				x += w*samples.x[ssize-i-1];
				y += w*samples.y[ssize-i-1];
				p += w*samples.pupil[ssize-i-1];
			}
			if (len%2==1) { 
				x += samples.x[start_fix+mid];
				y += samples.y[start_fix+mid];
				p += samples.pupil[start_fix+mid];
				++mid;
			}
			x /= mid+1;
			y /= mid+1;
			p /= mid+1;

			fix.start.push_back(start);
			fix.end.push_back(stop);
			fix.x.push_back(x);
			fix.y.push_back(y);
			fix.pupil.push_back(p);
		}

		/*
		 * Process one sample, respecting current fixations
		 */
		inline void sample(int time, float x, float y, float p){
			samples.time.push_back(time);
			samples.x.push_back(x);
			samples.y.push_back(y);
			samples.pupil.push_back(p);
		}
		#ifdef MEX
		/*
		 * Serialize to mxArray
		 */
		mxArray * asstruct(UINT32 start){
			if (this->samples.x.size() == 0)
				return mxCreateLogicalScalar( false );
			int dims[1] = { 1 };
			mxArray * result = mxCreateStructArray( 1, dims, EYE_COUNTF, eye_fields );
			mxSetFieldByNumber( result, 0, EYE_SACCADE, sac.toArray(start) );
			mxSetFieldByNumber( result, 0, EYE_FIXATION, fix.toArray(start) );
			mxSetFieldByNumber( result, 0, EYE_SAMPLES, samples.toArray(start) );
			mxSetFieldByNumber( result, 0, EYE_BLINK, blink.toArray(start) );
			mxSetFieldByNumber( result, 0, EYE_DRIFT, drift.toArray() );
			return result;
		}
		#endif
}; 

/*
 * Vector of Eyes. Must be init'ed to length=2!
 */
typedef vector<Eye> Eyes;
