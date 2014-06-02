
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "moses/Timer.h"
#include "moses/Hypothesis.h"
#include "moses/Manager.h"
#include "moses/StaticData.h"
#include "moses/Util.h"

#ifdef HAVE_PROTOBUF
#include "hypergraph.pb.h"
#endif

#include "TranslationTask.h"

namespace MosesCmd
{


    /** Enforce rounding */
    void fix(std::ostream& stream, size_t size)
    {
    stream.setf(std::ios::fixed);
    stream.precision(size);
    }

  /** Translate one sentence
   * gets called by main function implemented at end of this source file */
  void TranslationTask::Run() {
    // shorthand for "global data"
    const StaticData &staticData = StaticData::Instance();

    // input sentence
    Sentence sentence;

    // report wall time spent on translation
    Timer translationTime;
    translationTime.start();

    // report thread number
#if defined(WITH_THREADS) && defined(BOOST_HAS_PTHREADS)
    TRACE_ERR("Translating line " << m_lineNumber << "  in thread id " << pthread_self() << std::endl);
#endif


    // execute the translation
    // note: this executes the search, resulting in a search graph
    //       we still need to apply the decision rule (MAP, MBR, ...)
    Timer initTime;
    initTime.start();
    Manager manager(m_lineNumber, *m_source,staticData.GetSearchAlgorithm());
    VERBOSE(1, "Line " << m_lineNumber << ": Initialize search took " << initTime << " seconds total" << endl);
    manager.ProcessSentence();

    // we are done with search, let's look what we got
    Timer additionalReportingTime;
    additionalReportingTime.start();

    // output word graph
    if (m_wordGraphCollector) {
      ostringstream out;
      fix(out,PRECISION);
      manager.GetWordGraph(m_lineNumber, out);
      m_wordGraphCollector->Write(m_lineNumber, out.str());
    }

    // output search graph
    if (m_searchGraphCollector) {
      ostringstream out;
      fix(out,PRECISION);
      manager.OutputSearchGraph(m_lineNumber, out);
      m_searchGraphCollector->Write(m_lineNumber, out.str());

#ifdef HAVE_PROTOBUF
      if (staticData.GetOutputSearchGraphPB()) {
        ostringstream sfn;
        sfn << staticData.GetParam("output-search-graph-pb")[0] << '/' << m_lineNumber << ".pb" << ends;
        string fn = sfn.str();
        VERBOSE(2, "Writing search graph to " << fn << endl);
        fstream output(fn.c_str(), ios::trunc | ios::binary | ios::out);
        manager.SerializeSearchGraphPB(m_lineNumber, output);
      }
#endif
    }

    // Output search graph in HTK standard lattice format (SLF)
    if (m_outputSearchGraphSLF) {
      stringstream fileName;
      fileName << staticData.GetParam("output-search-graph-slf")[0] << "/" << m_lineNumber << ".slf";
      std::ofstream *file = new std::ofstream;
      file->open(fileName.str().c_str());
      if (file->is_open() && file->good()) {
        ostringstream out;
        fix(out,PRECISION);
        manager.OutputSearchGraphAsSLF(m_lineNumber, out);
        *file << out.str();
        file -> flush();
      } else {
        TRACE_ERR("Cannot output HTK standard lattice for line " << m_lineNumber << " because the output file is not open or not ready for writing" << std::endl);
      }
      delete file;
    }

    // Output search graph in hypergraph format for Kenneth Heafield's lazy hypergraph decoder
    if (m_outputSearchGraphHypergraph) {

      vector<string> hypergraphParameters = staticData.GetParam("output-search-graph-hypergraph");

      bool appendSuffix;
      if (hypergraphParameters.size() > 0 && hypergraphParameters[0] == "true") {
        appendSuffix = true;
      } else {
        appendSuffix = false;
      }

      string compression;
      if (hypergraphParameters.size() > 1) {
        compression = hypergraphParameters[1];
      } else {
        compression = "txt";
      }

      string hypergraphDir;
      if ( hypergraphParameters.size() > 2 ) {
        hypergraphDir = hypergraphParameters[2];
      } else {
        string nbestFile = staticData.GetNBestFilePath();
        if ( ! nbestFile.empty() && nbestFile!="-" && !boost::starts_with(nbestFile,"/dev/stdout") ) {
          boost::filesystem::path nbestPath(nbestFile);

          // In the Boost filesystem API version 2,
          //   which was the default prior to Boost 1.46,
          //   the filename() method returned a string.
          //
          // In the Boost filesystem API version 3,
          //   which is the default starting with Boost 1.46,
          //   the filename() method returns a path object.
          //
          // To get a string from the path object,
          //   the native() method must be called.
          //	  hypergraphDir = nbestPath.parent_path().filename()
          //#if BOOST_VERSION >= 104600
          //	    .native()
          //#endif
          //;

          // Hopefully the following compiles under all versions of Boost.
          //
          // If this line gives you compile errors,
          //   contact Lane Schwartz on the Moses mailing list
          hypergraphDir = nbestPath.parent_path().string();

        } else {
          stringstream hypergraphDirName;
          hypergraphDirName << boost::filesystem::current_path().string() << "/hypergraph";
          hypergraphDir = hypergraphDirName.str();
        }
      }

      if ( ! boost::filesystem::exists(hypergraphDir) ) {
        boost::filesystem::create_directory(hypergraphDir);
      }

      if ( ! boost::filesystem::exists(hypergraphDir) ) {
        TRACE_ERR("Cannot output hypergraphs to " << hypergraphDir << " because the directory does not exist" << std::endl);
      } else if ( ! boost::filesystem::is_directory(hypergraphDir) ) {
        TRACE_ERR("Cannot output hypergraphs to " << hypergraphDir << " because that path exists, but is not a directory" << std::endl);
      } else {
        stringstream fileName;
        fileName << hypergraphDir << "/" << m_lineNumber;
        if ( appendSuffix ) {
          fileName << "." << compression;
        }
        boost::iostreams::filtering_ostream *file 
	  = new boost::iostreams::filtering_ostream;

        if ( compression == "gz" ) {
          file->push( boost::iostreams::gzip_compressor() );
        } else if ( compression == "bz2" ) {
          file->push( boost::iostreams::bzip2_compressor() );
        } else if ( compression != "txt" ) {
          TRACE_ERR("Unrecognized hypergraph compression format (" 
		    << compression 
		    << ") - using uncompressed plain txt" << std::endl);
          compression = "txt";
        }

        file->push( boost::iostreams::file_sink(fileName.str(), ios_base::out) );

        if (file->is_complete() && file->good()) {
          fix(*file,PRECISION);
          manager.OutputSearchGraphAsHypergraph(m_lineNumber, *file);
          file -> flush();
        } else {
          TRACE_ERR("Cannot output hypergraph for line " << m_lineNumber 
		    << " because the output file " << fileName.str() 
		    << " is not open or not ready for writing" 
		    << std::endl);
        }
        file -> pop();
        delete file;
      }
    }
    additionalReportingTime.stop();

    // apply decision rule and output best translation(s)
    if (m_outputCollector) {
      ostringstream out;
      ostringstream debug;
      fix(debug,PRECISION);

      // all derivations - send them to debug stream
      if (staticData.PrintAllDerivations()) {
        additionalReportingTime.start();
        manager.PrintAllDerivations(m_lineNumber, debug);
        additionalReportingTime.stop();
      }

      Timer decisionRuleTime;
      decisionRuleTime.start();

      // MAP decoding: best hypothesis
      const Hypothesis* bestHypo = NULL;
      if (!staticData.UseMBR()) {
        bestHypo = manager.GetBestHypothesis();
        if (bestHypo) {
          if (StaticData::Instance().GetOutputHypoScore()) {
            out << bestHypo->GetTotalScore() << ' ';
          }
          if (staticData.IsPathRecoveryEnabled()) {
            OutputInput(out, bestHypo);
            out << "||| ";
          }
          if (staticData.GetParam("print-id").size() && Scan<bool>(staticData.GetParam("print-id")[0]) ) {
            out << m_source->GetTranslationId() << " ";
          }

	  if (staticData.GetReportSegmentation() == 2) {
	    manager.GetOutputLanguageModelOrder(out, bestHypo);
	  }
          OutputBestSurface(
            out,
            bestHypo,
            staticData.GetOutputFactorOrder(),
            staticData.GetReportSegmentation(),
            staticData.GetReportAllFactors());
          if (staticData.PrintAlignmentInfo()) {
            out << "||| ";
            OutputAlignment(out, bestHypo);
          }

          OutputAlignment(m_alignmentInfoCollector, m_lineNumber, bestHypo);
          IFVERBOSE(1) {
            debug << "BEST TRANSLATION: " << *bestHypo << endl;
          }
        } else {
          VERBOSE(1, "NO BEST TRANSLATION" << endl);
        }

        out << endl;
      }

      // MBR decoding (n-best MBR, lattice MBR, consensus)
      else {
        // we first need the n-best translations
        size_t nBestSize = staticData.GetMBRSize();
        if (nBestSize <= 0) {
          cerr << "ERROR: negative size for number of MBR candidate translations not allowed (option mbr-size)" << endl;
          exit(1);
        }
        TrellisPathList nBestList;
        manager.CalcNBest(nBestSize, nBestList,true);
        VERBOSE(2,"size of n-best: " << nBestList.GetSize() << " (" << nBestSize << ")" << endl);
        IFVERBOSE(2) {
          PrintUserTime("calculated n-best list for (L)MBR decoding");
        }

        // lattice MBR
        if (staticData.UseLatticeMBR()) {
          if (m_nbestCollector) {
            //lattice mbr nbest
            vector<LatticeMBRSolution> solutions;
            size_t n  = min(nBestSize, staticData.GetNBestSize());
            getLatticeMBRNBest(manager,nBestList,solutions,n);
            ostringstream out;
            OutputLatticeMBRNBest(out, solutions,m_lineNumber);
            m_nbestCollector->Write(m_lineNumber, out.str());
          } else {
            //Lattice MBR decoding
            vector<Word> mbrBestHypo = doLatticeMBR(manager,nBestList);
            OutputBestHypo(mbrBestHypo, m_lineNumber, staticData.GetReportSegmentation(),
                           staticData.GetReportAllFactors(),out);
            IFVERBOSE(2) {
              PrintUserTime("finished Lattice MBR decoding");
            }
          }
        }

        // consensus decoding
        else if (staticData.UseConsensusDecoding()) {
          const TrellisPath &conBestHypo = doConsensusDecoding(manager,nBestList);
          OutputBestHypo(conBestHypo, m_lineNumber,
                         staticData.GetReportSegmentation(),
                         staticData.GetReportAllFactors(),out);
          OutputAlignment(m_alignmentInfoCollector, m_lineNumber, conBestHypo);
          IFVERBOSE(2) {
            PrintUserTime("finished Consensus decoding");
          }
        }

        // n-best MBR decoding
        else {
          const Moses::TrellisPath &mbrBestHypo = doMBR(nBestList);
          OutputBestHypo(mbrBestHypo, m_lineNumber,
                         staticData.GetReportSegmentation(),
                         staticData.GetReportAllFactors(),out);
          OutputAlignment(m_alignmentInfoCollector, m_lineNumber, mbrBestHypo);
          IFVERBOSE(2) {
            PrintUserTime("finished MBR decoding");
          }
        }
      }

      // report best translation to output collector
      m_outputCollector->Write(m_lineNumber,out.str(),debug.str());

      decisionRuleTime.stop();
      VERBOSE(1, "Line " << m_lineNumber << ": Decision rule took " << decisionRuleTime << " seconds total" << endl);
    }

    additionalReportingTime.start();

    // output n-best list
    if (m_nbestCollector && !staticData.UseLatticeMBR()) {
      TrellisPathList nBestList;
      ostringstream out;
      manager.CalcNBest(staticData.GetNBestSize(), nBestList,staticData.GetDistinctNBest());
      OutputNBest(out, nBestList, staticData.GetOutputFactorOrder(), m_lineNumber,
                  staticData.GetReportSegmentation());
      m_nbestCollector->Write(m_lineNumber, out.str());
    }

    //lattice samples
    if (m_latticeSamplesCollector) {
      TrellisPathList latticeSamples;
      ostringstream out;
      manager.CalcLatticeSamples(staticData.GetLatticeSamplesSize(), latticeSamples);
      OutputNBest(out,latticeSamples, staticData.GetOutputFactorOrder(), m_lineNumber,
                  staticData.GetReportSegmentation());
      m_latticeSamplesCollector->Write(m_lineNumber, out.str());
    }

    // detailed translation reporting
    if (m_detailedTranslationCollector) {
      ostringstream out;
      fix(out,PRECISION);
      TranslationAnalysis::PrintTranslationAnalysis(out, manager.GetBestHypothesis());
      m_detailedTranslationCollector->Write(m_lineNumber,out.str());
    }

    //list of unknown words
    if (m_unknownsCollector) {
      const vector<const Phrase*>& unknowns = manager.getSntTranslationOptions()->GetUnknownSources();
      ostringstream out;
      for (size_t i = 0; i < unknowns.size(); ++i) {
        out << *(unknowns[i]);
      }
      out << endl;
      m_unknownsCollector->Write(m_lineNumber, out.str());
    }

    // report additional statistics
    manager.CalcDecoderStatistics();
    VERBOSE(1, "Line " << m_lineNumber << ": Additional reporting took " << additionalReportingTime << " seconds total" << endl);
    VERBOSE(1, "Line " << m_lineNumber << ": Translation took " << translationTime << " seconds total" << endl);
    IFVERBOSE(2) {
      PrintUserTime("Sentence Decoding Time:");
    }
  }

  TranslationTask::~TranslationTask() {
    delete m_source;
  }

} //end namespace
