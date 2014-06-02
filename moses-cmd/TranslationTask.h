// $Id$
//
//

#include <exception>
#include <fstream>
#include <sstream>
#include <vector>

#include "util/usage.hh"

#include "TranslationAnalysis.h"
#include "IOWrapper.h"
#include "mbr.h"

#include "moses/ThreadPool.h"
#include "moses/OutputCollector.h"

using namespace std;
using namespace Moses;
using namespace MosesCmd;

namespace MosesCmd
{

// output floats with five significant digits
static const size_t PRECISION = 3;
/** Enforce rounding */
void fix(std::ostream& stream, size_t size);

/** Translates a sentence.
  * - calls the search (Manager)
  * - applies the decision rule
  * - outputs best translation and additional reporting
  **/
class TranslationTask : public Task
{

public:

  TranslationTask(size_t lineNumber,
                  InputType* source, OutputCollector* outputCollector, OutputCollector* nbestCollector,
                  OutputCollector* latticeSamplesCollector,
                  OutputCollector* wordGraphCollector, OutputCollector* searchGraphCollector,
                  OutputCollector* detailedTranslationCollector,
                  OutputCollector* alignmentInfoCollector,
                  OutputCollector* unknownsCollector,
                  bool outputSearchGraphSLF,
                  bool outputSearchGraphHypergraph) :
    m_source(source), m_lineNumber(lineNumber),
    m_outputCollector(outputCollector), m_nbestCollector(nbestCollector),
    m_latticeSamplesCollector(latticeSamplesCollector),
    m_wordGraphCollector(wordGraphCollector), m_searchGraphCollector(searchGraphCollector),
    m_detailedTranslationCollector(detailedTranslationCollector),
    m_alignmentInfoCollector(alignmentInfoCollector),
    m_unknownsCollector(unknownsCollector),
    m_outputSearchGraphSLF(outputSearchGraphSLF),
    m_outputSearchGraphHypergraph(outputSearchGraphHypergraph) {}

  /** Translate one sentence
   * gets called by main function implemented at end of this Main.cpp */
  void Run();

  ~TranslationTask();

private:
  InputType* m_source;
  size_t m_lineNumber;
  OutputCollector* m_outputCollector;
  OutputCollector* m_nbestCollector;
  OutputCollector* m_latticeSamplesCollector;
  OutputCollector* m_wordGraphCollector;
  OutputCollector* m_searchGraphCollector;
  OutputCollector* m_detailedTranslationCollector;
  OutputCollector* m_alignmentInfoCollector;
  OutputCollector* m_unknownsCollector;
  bool m_outputSearchGraphSLF;
  bool m_outputSearchGraphHypergraph;
  std::ofstream *m_alignmentStream;


};

} //end namespace
