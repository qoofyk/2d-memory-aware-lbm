#include <papi.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_EVENTS 2
#define ERROR_RETURN(retval) { fprintf(stderr, "Error %d %s:line %d: \n", retval,__FILE__,__LINE__);  exit(retval); }

int main(int argc, char *argv[]){
	int EventSet = PAPI_NULL;
	int native = 0x0;
	char error_str[PAPI_MAX_STR_LEN];
	long long values[2];
	int number, retval;

	/* Initialize the PAPI library */
	if((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT )
      ERROR_RETURN(retval);
     

	/* Creating the eventset */              
   if ( (retval = PAPI_create_eventset(&EventSet)) != PAPI_OK)
      ERROR_RETURN(retval);

	/* Add Total Cycles event to the EventSet */
   if ( (retval = PAPI_add_event(EventSet, PAPI_TOT_CYC)) != PAPI_OK)
      ERROR_RETURN(retval);

	/* Add Total Instructions Executed to our EventSet */
	if ((retval = PAPI_add_event(EventSet, PAPI_TOT_INS)) != PAPI_OK) {
		PAPI_perror(error_str);
		fprintf(stderr,"PAPI_add_event %d: %s\n", retval, error_str);
		exit(1);
	}

	/* get the number of events in the event set */
   number = 0;
   if ( (retval = PAPI_list_events(EventSet, NULL, &number)) != PAPI_OK)
      ERROR_RETURN(retval);

   printf("There are %d events in the event set\n", number);

	// /* Add LLC_REFERENCES native event */
	// retval = PAPI_event_name_to_code("LLC_REFERENCES", &native);
	// if ((retval = PAPI_add_event(&EventSet, native)) != PAPI_OK) {
	// 	fprintf(stderr, "PAPI_add_event error %d: %s\n", retval, PAPI_strerror(retval));
	// 	exit(1);
	// }

	/* Start counting */
	if ((retval = PAPI_start(EventSet)) != PAPI_OK){
		fprintf(stderr, "PAPI_start error %d: %s\n", retval, PAPI_strerror(retval));
		exit(1);
	}

	/* Stop the counting of events in the Event Set */
	if ((retval = PAPI_stop(EventSet, values)) != PAPI_OK){
		fprintf(stderr, "PAPI_stop error %d: %s\n", retval, PAPI_strerror(retval));
		exit(1);
	}

	/* Remove event: We are going to take the PAPI_TOT_INS from the eventset */
   if( (retval = PAPI_remove_event(EventSet, PAPI_TOT_INS)) != PAPI_OK)
      ERROR_RETURN(retval);
   printf("Removing PAPI_TOT_INS from the eventset\n"); 

   /* Now we list how many events are left on the event set */
   number = 0;
   if ((retval=PAPI_list_events(EventSet, NULL, &number))!= PAPI_OK)
      ERROR_RETURN(retval);

   printf("There is only %d event left in the eventset now\n", number);

   /* free the resources used by PAPI */
   PAPI_shutdown();

	return 1;
}