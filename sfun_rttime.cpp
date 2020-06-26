#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME  sfun_rttime

#define TIME_SCALE_FACTOR(S) ssGetSFcnParam(S,0)

/* Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions. */
#include <simstruc.h>

#if defined(_WIN32)
/* Include the windows SDK header for handling time functions. */
#include <windows.h>
#include <math.h>

struct TimerData
{
	TimerData() {}
	_LARGE_INTEGER frequency;
	_LARGE_INTEGER previous_time;
};

class Timer
{
public:
	Timer() : m_data(NULL)
	{
		Init();
	}

	~Timer()
	{
		delete m_data;
	}
	double Event()
	{
		_LARGE_INTEGER old_time = m_data->previous_time;
		QueryPerformanceCounter(&m_data->previous_time);
		double DeltaT = double(m_data->previous_time.QuadPart-old_time.QuadPart)/double(m_data->frequency.QuadPart);
		return DeltaT;
	}
	double TimeSinceLastEvent()
	{
		_LARGE_INTEGER current_time;
		QueryPerformanceCounter(&current_time);
		return double(current_time.QuadPart-m_data->previous_time.QuadPart)/double(m_data->frequency.QuadPart);
	}
	void Init()
	{
		delete m_data;
		m_data = new TimerData();
		QueryPerformanceFrequency(&m_data->frequency);
		QueryPerformanceCounter(&m_data->previous_time);
	}
private:
	TimerData* m_data;
};

void SleepSeconds(double t)
{
	Sleep(unsigned long(t*1000.0));
}

static Timer masterTimer;
static int stepNumber = 0;
static double initSimulinkTime = 0;

#else
/* Include the standard ANSI C header for handling time functions. */
#include <time.h>

/* Function of the high performance counter (in seconds). */
__inline double hightimer()
{    
    return (double)clock()/CLOCKS_PER_SEC;
}	/* end hightimer */
#endif

static void mdlInitializeSizes(SimStruct *S)
{
   ssSetNumSFcnParams(S, 1);  /* Number of expected parameters */
   if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return;
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 1);
   if (!ssSetNumInputPorts(S, 0)) return;
   if (!ssSetNumOutputPorts(S, 1)) return;
   ssSetOutputPortWidth(S, 0, 1);
   ssSetNumSampleTimes(S, 1);
   ssSetNumRWork(S, 1);
   ssSetNumIWork(S, 0);
   ssSetNumPWork(S, 0);
   ssSetNumModes(S, 0);
   ssSetNumNonsampledZCs(S, 0);
   ssSetOptions(S, 0);
}

#define MDL_INITIALIZE_SAMPLE_TIMES
static void mdlInitializeSampleTimes(SimStruct *S)
{
   ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
   ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START
static void mdlStart(SimStruct *S)
{
   ssSetRWorkValue(S,0,ssGetTStart(S));

   masterTimer.Event();
   stepNumber = 0;
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
	double currentSimTime = ssGetT(S);
	   //std::cout << currentSimTime << std::endl;
	   //mexPrintf("time: %lf\n", currentSimTime);
	   
	   
   if (stepNumber++<10)
   {
		masterTimer.Event();
		initSimulinkTime = currentSimTime;
   }
   else
   {
  
	   const double *scaleFactor = mxGetPr(TIME_SCALE_FACTOR(S));
	   while (((masterTimer.TimeSinceLastEvent()-initSimulinkTime)*(*scaleFactor))<currentSimTime)
			int toto = 0;
   }
   /*while (masterTimer.TimeSinceLastEvent()<0.002)
	   int toto = 1;
   masterTimer.Event();*/
}

static void mdlTerminate(SimStruct *S)
{
    UNUSED_ARG(S); /* unused input argument */
}

/* Required S-function trailer */
#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
