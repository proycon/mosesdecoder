from libcpp.string cimport string
from libcpp.vector cimport vector

class MosesError(Exception):
    pass

cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        pass

    extern ostream cout
    extern ostream cerr

cdef extern from 'OutputCollector.h' namespace 'Moses':
    cdef cppclass OutputCollector:
        OutputCollector(ostream*, ostream*)

cdef extern from "InputType.h" namespace 'Moses':
    cdef cppclass InputType:
        pass

cdef extern from "Parameter.h" namespace 'Moses':
    cdef cppclass Parameter:
        bool LoadParam(int argc, char** argv)

cdef extern from 'TranslationTask.h' namespace 'MosesCmd':
    cdef cppclass TranslationTask(int lineNumber,
                  InputType* source, OutputCollector* outputCollector, OutputCollector* nbestCollector,
                  OutputCollector* latticeSamplesCollector,
                  OutputCollector* wordGraphCollector, OutputCollector* searchGraphCollector,
                  OutputCollector* detailedTranslationCollector,
                  OutputCollector* alignmentInfoCollector,
                  OutputCollector* unknownsCollector,
                  bool outputSearchGraphSLF,
                  bool outputSearchGraphHypergraph):
        void Run() nogil +MosesError

    fix(ostream*, int)
