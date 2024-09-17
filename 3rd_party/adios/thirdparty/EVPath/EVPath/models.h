#include <cod.h>
#include <evpath.h>


// fake enum types
// TODO: replace with real enums at a future date when it is implemented in EVPath
extern const int MODEL_NO_FIT;
extern const int MODEL_FFT;

// TODO: this is probably a good struct to hold a model, but err_term should probably be a pointer to allow for multiple errors
// holds the information that describe the parameters of any model of any data
typedef struct _model {
	int type;
	int num_terms;
	double* terms;
	double err_term;
} model, *model_ptr;

// function prototypes
int EVmodel_init(CManager cm);