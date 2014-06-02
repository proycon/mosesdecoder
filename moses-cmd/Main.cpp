// $Id: MainMT.cpp 3045 2010-04-05 13:07:29Z hieuhoang1972 $

/***********************************************************************
Moses - factored phrase-based language decoder
Copyright (C) 2009 University of Edinburgh

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
***********************************************************************/

/**
 * Moses main, for single-threaded and multi-threaded.
 **/

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <exception>
#include <fstream>
#include <sstream>
#include <vector>

#include "util/usage.hh"

#ifdef WIN32
// Include Visual Leak Detector
//#include <vld.h>
#endif

#include "TranslationAnalysis.h"
#include "IOWrapper.h"
#include "mbr.h"


#include "moses/StaticData.h"
#include "moses/Util.h"
#include "moses/ThreadPool.h"
#include "moses/TranslationModel/PhraseDictionary.h"
#include "moses/FF/StatefulFeatureFunction.h"
#include "moses/FF/StatelessFeatureFunction.h"

#ifdef HAVE_PROTOBUF
#include "hypergraph.pb.h"
#endif

//Contains the main Run() method
#include "TranslationTask.h" 

using namespace std;
using namespace Moses;
using namespace MosesCmd;

namespace MosesCmd
{


static void PrintFeatureWeight(const FeatureFunction* ff)
{
  cout << ff->GetScoreProducerDescription() << "=";
  size_t numScoreComps = ff->GetNumScoreComponents();
  vector<float> values = StaticData::Instance().GetAllWeights().GetScoresForProducer(ff);
  for (size_t i = 0; i < numScoreComps; ++i) {
    cout << " " << values[i];
  }
  cout << endl;
}

static void ShowWeights()
{
  //TODO: Find a way of ensuring this order is synced with the nbest
  fix(cout,6);
  const vector<const StatelessFeatureFunction*>& slf = StatelessFeatureFunction::GetStatelessFeatureFunctions();
  const vector<const StatefulFeatureFunction*>& sff = StatefulFeatureFunction::GetStatefulFeatureFunctions();

  for (size_t i = 0; i < sff.size(); ++i) {
    const StatefulFeatureFunction *ff = sff[i];
    if (ff->IsTuneable()) {
      PrintFeatureWeight(ff);
    }
    else {
      cout << ff->GetScoreProducerDescription() << " UNTUNEABLE" << endl;
    }
  }
  for (size_t i = 0; i < slf.size(); ++i) {
    const StatelessFeatureFunction *ff = slf[i];
    if (ff->IsTuneable()) {
      PrintFeatureWeight(ff);
    }
    else {
      cout << ff->GetScoreProducerDescription() << " UNTUNEABLE" << endl;
    }
  }
}

size_t OutputFeatureWeightsForHypergraph(size_t index, const FeatureFunction* ff, std::ostream &outputSearchGraphStream)
{
  size_t numScoreComps = ff->GetNumScoreComponents();
  if (numScoreComps != 0) {
    vector<float> values = StaticData::Instance().GetAllWeights().GetScoresForProducer(ff);
    if (numScoreComps > 1) {
      for (size_t i = 0; i < numScoreComps; ++i) {
        outputSearchGraphStream << ff->GetScoreProducerDescription()
                                << i
                                << "=" << values[i] << endl;
      }
    } else {
      outputSearchGraphStream << ff->GetScoreProducerDescription()
                              << "=" << values[0] << endl;
    }
    return index+numScoreComps;
  } else {
    UTIL_THROW2("Sparse features are not yet supported when outputting hypergraph format");
  }
}

void OutputFeatureWeightsForHypergraph(std::ostream &outputSearchGraphStream)
{
  outputSearchGraphStream.setf(std::ios::fixed);
  outputSearchGraphStream.precision(6);

  const vector<const StatelessFeatureFunction*>& slf =StatelessFeatureFunction::GetStatelessFeatureFunctions();
  const vector<const StatefulFeatureFunction*>& sff = StatefulFeatureFunction::GetStatefulFeatureFunctions();
  size_t featureIndex = 1;
  for (size_t i = 0; i < sff.size(); ++i) {
    featureIndex = OutputFeatureWeightsForHypergraph(featureIndex, sff[i], outputSearchGraphStream);
  }
  for (size_t i = 0; i < slf.size(); ++i) {
    /*
    if (slf[i]->GetScoreProducerWeightShortName() != "u" &&
          slf[i]->GetScoreProducerWeightShortName() != "tm" &&
          slf[i]->GetScoreProducerWeightShortName() != "I" &&
          slf[i]->GetScoreProducerWeightShortName() != "g")
    */
    {
      featureIndex = OutputFeatureWeightsForHypergraph(featureIndex, slf[i], outputSearchGraphStream);
    }
  }
  const vector<PhraseDictionary*>& pds = PhraseDictionary::GetColl();
  for( size_t i=0; i<pds.size(); i++ ) {
    featureIndex = OutputFeatureWeightsForHypergraph(featureIndex, pds[i], outputSearchGraphStream);
  }
  const vector<GenerationDictionary*>& gds = GenerationDictionary::GetColl();
  for( size_t i=0; i<gds.size(); i++ ) {
    featureIndex = OutputFeatureWeightsForHypergraph(featureIndex, gds[i], outputSearchGraphStream);
  }

}


} //namespace


/** main function of the command line version of the decoder **/
int main(int argc, char** argv)
{
  try {

#ifdef HAVE_PROTOBUF
    GOOGLE_PROTOBUF_VERIFY_VERSION;
#endif
    
    // echo command line, if verbose
    IFVERBOSE(1) {
      TRACE_ERR("command: ");
      for(int i=0; i<argc; ++i) TRACE_ERR(argv[i]<<" ");
      TRACE_ERR(endl);
    }

    // set number of significant decimals in output
    fix(cout,PRECISION);
    fix(cerr,PRECISION);

    // load all the settings into the Parameter class
    // (stores them as strings, or array of strings)
    Parameter params;
    if (!params.LoadParam(argc,argv)) {
      exit(1);
    }


    // initialize all "global" variables, which are stored in StaticData
    // note: this also loads models such as the language model, etc.
    if (!StaticData::LoadDataStatic(&params, argv[0])) {
      exit(1);
    }

    // setting "-show-weights" -> just dump out weights and exit
    if (params.isParamSpecified("show-weights")) {
      ShowWeights();
      exit(0);
    }

    // shorthand for accessing information in StaticData
    const StaticData& staticData = StaticData::Instance();


    //initialise random numbers
    srand(time(NULL));

    // set up read/writing class
    IOWrapper* ioWrapper = GetIOWrapper(staticData);
    if (!ioWrapper) {
      cerr << "Error; Failed to create IO object" << endl;
      exit(1);
    }

    // check on weights
    const ScoreComponentCollection& weights = staticData.GetAllWeights();
    IFVERBOSE(2) {
      TRACE_ERR("The global weight vector looks like this: ");
      TRACE_ERR(weights);
      TRACE_ERR("\n");
    }
    if (staticData.GetOutputSearchGraphHypergraph()) {
      ofstream* weightsOut = new std::ofstream;
      stringstream weightsFilename;
      if (staticData.GetParam("output-search-graph-hypergraph").size() > 3) {
        weightsFilename << staticData.GetParam("output-search-graph-hypergraph")[3];
      } else {
        string nbestFile = staticData.GetNBestFilePath();
        if ( ! nbestFile.empty() && nbestFile!="-" && !boost::starts_with(nbestFile,"/dev/stdout") ) {
          boost::filesystem::path nbestPath(nbestFile);
          weightsFilename << nbestPath.parent_path().filename() << "/weights";
        } else {
          weightsFilename << boost::filesystem::current_path().string() << "/hypergraph/weights";
        }
      }
      boost::filesystem::path weightsFilePath(weightsFilename.str());
      if ( ! boost::filesystem::exists(weightsFilePath.parent_path()) ) {
        boost::filesystem::create_directory(weightsFilePath.parent_path());
      }
      TRACE_ERR("The weights file is " << weightsFilename.str() << "\n");
      weightsOut->open(weightsFilename.str().c_str());
      OutputFeatureWeightsForHypergraph(*weightsOut);
      weightsOut->flush();
      weightsOut->close();
      delete weightsOut;
    }


    // initialize output streams
    // note: we can't just write to STDOUT or files
    // because multithreading may return sentences in shuffled order
    auto_ptr<OutputCollector> outputCollector; // for translations
    auto_ptr<OutputCollector> nbestCollector;  // for n-best lists
    auto_ptr<OutputCollector> latticeSamplesCollector; //for lattice samples
    auto_ptr<ofstream> nbestOut;
    auto_ptr<ofstream> latticeSamplesOut;
    size_t nbestSize = staticData.GetNBestSize();
    string nbestFile = staticData.GetNBestFilePath();
    bool output1best = true;
    if (nbestSize) {
      if (nbestFile == "-" || nbestFile == "/dev/stdout") {
        // nbest to stdout, no 1-best
        nbestCollector.reset(new OutputCollector());
        output1best = false;
      } else {
        // nbest to file, 1-best to stdout
        nbestOut.reset(new ofstream(nbestFile.c_str()));
        if (!nbestOut->good()) {
          TRACE_ERR("ERROR: Failed to open " << nbestFile << " for nbest lists" << endl);
          exit(1);
        }
        nbestCollector.reset(new OutputCollector(nbestOut.get()));
      }
    }
    size_t latticeSamplesSize = staticData.GetLatticeSamplesSize();
    string latticeSamplesFile = staticData.GetLatticeSamplesFilePath();
    if (latticeSamplesSize) {
      if (latticeSamplesFile == "-" || latticeSamplesFile == "/dev/stdout") {
        latticeSamplesCollector.reset(new OutputCollector());
        output1best = false;
      } else {
        latticeSamplesOut.reset(new ofstream(latticeSamplesFile.c_str()));
        if (!latticeSamplesOut->good()) {
          TRACE_ERR("ERROR: Failed to open " << latticeSamplesFile << " for lattice samples" << endl);
          exit(1);
        }
        latticeSamplesCollector.reset(new OutputCollector(latticeSamplesOut.get()));
      }
    }
    if (output1best) {
      outputCollector.reset(new OutputCollector());
    }

    // initialize stream for word graph (aka: output lattice)
    auto_ptr<OutputCollector> wordGraphCollector;
    if (staticData.GetOutputWordGraph()) {
      wordGraphCollector.reset(new OutputCollector(&(ioWrapper->GetOutputWordGraphStream())));
    }

    // initialize stream for search graph
    // note: this is essentially the same as above, but in a different format
    auto_ptr<OutputCollector> searchGraphCollector;
    if (staticData.GetOutputSearchGraph()) {
      searchGraphCollector.reset(new OutputCollector(&(ioWrapper->GetOutputSearchGraphStream())));
    }

    // initialize stram for details about the decoder run
    auto_ptr<OutputCollector> detailedTranslationCollector;
    if (staticData.IsDetailedTranslationReportingEnabled()) {
      detailedTranslationCollector.reset(new OutputCollector(&(ioWrapper->GetDetailedTranslationReportingStream())));
    }

    // initialize stram for word alignment between input and output
    auto_ptr<OutputCollector> alignmentInfoCollector;
    if (!staticData.GetAlignmentOutputFile().empty()) {
      alignmentInfoCollector.reset(new OutputCollector(ioWrapper->GetAlignmentOutputStream()));
    }

    //initialise stream for unknown (oov) words
    auto_ptr<OutputCollector> unknownsCollector;
    auto_ptr<ofstream> unknownsStream;
    if (!staticData.GetOutputUnknownsFile().empty()) {
      unknownsStream.reset(new ofstream(staticData.GetOutputUnknownsFile().c_str()));
      if (!unknownsStream->good()) {
        TRACE_ERR("Unable to open " << staticData.GetOutputUnknownsFile() << " for unknowns");
        exit(1);
      }
      unknownsCollector.reset(new OutputCollector(unknownsStream.get()));
    }

#ifdef WITH_THREADS
    ThreadPool pool(staticData.ThreadCount());
#endif

    // main loop over set of input sentences
    InputType* source = NULL;
    size_t lineCount = staticData.GetStartTranslationId();
    while(ReadInput(*ioWrapper,staticData.GetInputType(),source)) {
      IFVERBOSE(1) {
        ResetUserTime();
      }
      // set up task of translating one sentence
      TranslationTask* task =
        new TranslationTask(lineCount,source, outputCollector.get(),
                            nbestCollector.get(),
                            latticeSamplesCollector.get(),
                            wordGraphCollector.get(),
                            searchGraphCollector.get(),
                            detailedTranslationCollector.get(),
                            alignmentInfoCollector.get(),
                            unknownsCollector.get(),
                            staticData.GetOutputSearchGraphSLF(),
                            staticData.GetOutputSearchGraphHypergraph());
      // execute task
#ifdef WITH_THREADS
      pool.Submit(task);
#else
      task->Run();
      delete task;
#endif

      source = NULL; //make sure it doesn't get deleted
      ++lineCount;
    }

    // we are done, finishing up
#ifdef WITH_THREADS
    pool.Stop(true); //flush remaining jobs
#endif

    delete ioWrapper;
    FeatureFunction::Destroy();

  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  IFVERBOSE(1) util::PrintUsage(std::cerr);

#ifndef EXIT_RETURN
  //This avoids that destructors are called (it can take a long time)
  exit(EXIT_SUCCESS);
#else
  return EXIT_SUCCESS;
#endif
}
