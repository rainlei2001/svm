// Disable the annoying debug warning on NT.
#ifdef WIN32
#pragma warning(disable:4786)
#endif

#include <math.h>
#include <vector>

// Applications headers
#include <apps/ExprArgs.h>
#include <apps/AppsUtil.h>
#include <apps/svm.h>

// RPAS headers
#include <rpas/SpecialExpression.h>
#include <rpas/Parser.h>
#include <rpas/Intersection.h>
#include <rpas/MeasureType.h>
#include <rpas/VarParseNode.h>
#include <rpas/ConstParseNode.h>
#include <rpas/MeasureWrapper.h>
#include <rpas/IllegalParse.h>
#include <rpas/ExpressionExceptions.h>
#include <rpas/Array.h>
#include <rpas/UPLogicalIterator.h>
#include <rpas/UPPopulatedIterator.h>
#include <rpas/UPPopVecIterator.h>
#include <rpas/ArrayMap.h>
#include <rpas/Dimension.h>
#include <rpas/Domain.h>
#include <rpas/SimpleArrayWrapper.h>
#include <rpas/Logger.h>
#include <rpas/ArrayExceptions.h>
#include <rdf/RDFSysConstants.h>
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

using namespace std;

namespace rdf
{
	RETEK_EXCEPTION_1(RdfSvmExprGeneralException, rpas::String)
	RETEK_EXCEPTION_1(RdfSvmExprException, rpas::String)

	const int _RdfSvmExprOutputArgsNum = 3;
	int _RdfSvmExprFixedInputNum = 0;
	int _dynamicFeaturesNum=0;

	int _totalArgNum=0;

	enum { PREDICT=1,OPT_PARAMETER,PRE_OPT };

   //----------------------------------------------------------------------
   //
   // RdfSvmExprArgs
   //
   //-----------------------------------------------------------------------

	class RdfSvmExprArgs : public apps::ExprArgs
	{
	public:
		// define input/output enum for RdfSvmExpr
		typedef enum
		{
		//output
         SVM_OUTPUT=0,
		 SVM_OPTGAMMA,
		 SVM_OPTCONST,

		//Input
		 RUN_MASK,    //mask to run
 		 SVM_TARGET,    //data value (y)
         SVM_TYPE,     //  C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR
		 SVM_OPERATION, //1-->preict,2-->optimize parameter ,3-->predict + opt parameters   
         SVM_KERNEL_TYPE,//{ LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED }; /* kernel_type *
         SVM_DEGREE, //* for poly */
		 SVM_GAMMA,  /* for poly/rbf/sigmoid */
		 SVM_COEF, /* for poly/sigmoid */
		 SVM_CACHESIZE, /* in MB */
		 SVM_EPS,/* stopping criteria */
		 SVM_CONST, /* for C_SVC, EPSILON_SVR and NU_SVR */
		 SVM_NU, /* for NU_SVC, ONE_CLASS, and NU_SVR */
		 SVM_EPSION,/* for EPSILON_SVR */
		 SVM_SHRINKING,/* use the shrinking heuristics */
		 SVM_NFOLDER, /*n folder for parameter optimizatoin */
		 SVM_ALPHA, /*the alpha value for optimization parameters */
		 SVM_TRAIN_BEGIN, /*begin index of the trainning window */
		 SVM_TRAIN_END, /*end index of the training window */
		 SVM_PREDICT_BEGIN, /*begin index of the predict hirizon */
		 SVM_PREDICT_END, /*end index of the predict hirizon*/
		 DATA_HIER,
		 TOTAL_ARGS
		} RdfSvmExprMeasureLabel;


		void fillMap()
		{
			//Output	
			//optional=false, input=false 
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_OUTPUT"), SVM_OUTPUT));
            _paraVec.push_back(apps::ExprParaInfo(false, false, apps::ExprParaInfo::MPI_MEASURE,
                                                  apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::real()));   
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_OPTGAMMA"), SVM_OPTGAMMA));
            _paraVec.push_back(apps::ExprParaInfo(false, false, apps::ExprParaInfo::MPI_MEASURE,
                                                  apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::real()));   
			
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_OPTCONST"), SVM_OPTCONST));
            _paraVec.push_back(apps::ExprParaInfo(false, false, apps::ExprParaInfo::MPI_MEASURE,
                                                  apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::real()));   
			


			//input
			//optional=false, input=true 
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("RUN_MASK"), RUN_MASK));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
                                                  apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::boolean()));   

			//SVM_TARGET
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_TARGET"), SVM_TARGET));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
                                                  apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::real())); 
			//SVM_TYPE
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_TYPE"), SVM_TYPE));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
                                                  apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::integer())); 
 
			//SVM_OPERATION
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_OPERATION"), SVM_OPERATION));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
                                                  apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::integer())); 
 


			//SVM_KERNEL_TYPE
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_KERNEL_TYPE"), SVM_KERNEL_TYPE));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
                                                  apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::integer()));  
						//SVM_DEGREE
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_DEGREE"), SVM_DEGREE));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
                                                  apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::integer()));  
						//SVM_GAMMA
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_GAMMA"), SVM_GAMMA));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
                                                  apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::real()));  
						//SVM_COEF
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_COEF"), SVM_COEF));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
                                                  apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::real()));  
						//SVM_CACHESIZE
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_CACHESIZE"), SVM_CACHESIZE));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
												 apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::integer()));  
			//SVM_EPS
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_EPS"), SVM_EPS));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
														apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::real()));  
						//SVM_CONST
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_CONST"), SVM_CONST));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
														apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::real()));  
						//SVM_NU
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_NU"), SVM_NU));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
														apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::real()));  
						//SVM_EPSION
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_EPSION"), SVM_EPSION));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
														apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::real()));  
						//SVM_SHRINKING
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_SHRINKING"), SVM_SHRINKING));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
														apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::integer())); 
				//SVM_NFOLDER
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_NFOLDER"), SVM_NFOLDER));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
														apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::integer())); 

				//SVM_ALPHA
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_ALPHA"), SVM_ALPHA));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
														apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::real())); 




									//SVM_TRAIN_BEGIN
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_TRAIN_BEGIN"), SVM_TRAIN_BEGIN));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
														apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::integer())); 
									//SVM_TRAIN_END
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_TRAIN_END"), SVM_TRAIN_END));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
														apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::integer())); 
									//SVM_PREDICT_BEGIN
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_PREDICT_BEGIN"), SVM_PREDICT_BEGIN));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
														apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::integer())); 
									//SVM_PREDICT_END
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_PREDICT_END"), SVM_PREDICT_END));
            _paraVec.push_back(apps::ExprParaInfo(false, true, apps::ExprParaInfo::MPI_MEASURE,
														apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::integer())); 
			 // DATA_HIER
			_labelMap.insert(Label2EnumMap::value_type(rpas::String("DATA_HIER"), DATA_HIER));
			 _paraVec.push_back(apps::ExprParaInfo(false, true, 
                                            apps::ExprParaInfo::MPI_SCALAR,
                                            apps::ExprParaInfo::CPI_CLND_DIM_PROHIBITED,
                                            rpas::MeasureType::string())); 




      }

      void fillSvmMap(int num) {
         fillMap();
         _dynamicFeaturesNum = num;
         for (int i=1 ; i<= num; i++ ) { 
            _labelMap.insert(Label2EnumMap::value_type(rpas::String("SVM_FEATURE"+rpas::String::intToString(i)),
               DATA_HIER+i));
            _paraVec.push_back(apps::ExprParaInfo(true, true, apps::ExprParaInfo::MPI_MEASURE,
               apps::ExprParaInfo::CPI_CLND_DIM_ALLOWED, rpas::MeasureType::real()));
         }
         _totalArgNum = _RdfSvmExprOutputArgsNum + _RdfSvmExprFixedInputNum + _dynamicFeaturesNum;
		 
      }

   };//class RdfSvmExprArgs

   
   
   
   //----------------------------------------------------------------------
   //
   // RdfSvmExpr
   //
   //-----------------------------------------------------------------------
   
   
	class RdfSvmExpr:public rpas::SpecialExpression
	{
      //Data
	private:
      RdfSvmExprArgs  _args;
	  rpas::ExpressionEvalMethod _ExpressionEvalMode;
  	   std::vector<rpas::Array> _InputArrays;
	 
      vector<rpas::SimpleArrayWrapper*> _OutputArrayVec;

      rpas::String _dataHier;
	  rpas::ArrayMap* _runMaskToTarget;

	  int _dataDimIndexInTarget ;
	  int _dataDimSize;
      int _feaureNaValue;
	  int _forecastLength;
	  svm_model * _model;
	  svm_problem _prob;		
      RdfSvmExpr();
      RdfSvmExpr(const RdfSvmExpr& src);            // Copy constructor
      RdfSvmExpr& operator=(const RdfSvmExpr& rhs); // Assignment operator

	   
	  //-----Function----------------------------------------------------
      //-----------------------------------------------------------------

 
public:

	  void cleanup()
      {

			for (int i = 0; i < _RdfSvmExprOutputArgsNum; i++)
			{
					delete _OutputArrayVec[i];
			}

			if (_runMaskToTarget !=0)
			{
				delete _runMaskToTarget;

			}
			
			if (_model !=0)
			{
				svm_free_and_destroy_model(&_model);
			}

			  for (int i=0;i<_prob.l;i++)
			  {
				  if (_prob.x[i] !=0)
				  {
					  free(_prob.x[i]);
				  }
			  }


			free(_prob.x);
			free(_prob.y);

      }

      virtual ~RdfSvmExpr() 
      {
      }

      RdfSvmExpr(const rpas::String& infixString,
         const rpas::String& expressionName, 
         const std::vector<rpas::ParseNode *>& lhs,
         const std::vector<rpas::ParseNode *>& rhs)
         : SpecialExpression(infixString, expressionName, lhs, rhs)
      {

			_runMaskToTarget = 0;
			_model=new svm_model();



      }

      //-----------validation------------------------------
      virtual void specialExpressionValidateHook()
      {
		  _RdfSvmExprFixedInputNum = RdfSvmExprArgs::TOTAL_ARGS - _RdfSvmExprOutputArgsNum;

         int featureNum = (_rhs.size()-_RdfSvmExprFixedInputNum);
         _args.fillSvmMap(featureNum);

         int reminder = (_rhs.size()-_RdfSvmExprFixedInputNum -featureNum );
         if (reminder > 0) {
            usage();
            throw rpas::IllegalParse("The number of inputs are not correct");
         }

         if (featureNum <1) {
            usage();
            throw rpas::IllegalParse("There is no feature inputs");
         }
         _args.validateParseNodesInSpecialExpression(_lhs);
         _args.validateParseNodesInSpecialExpression(_rhs);
         if (!_args.validInput())
         {
            usage();
            throw rpas::IllegalParse("Necessary inputs are missing");
         }
         if (!_args.validOutput())
         {
            usage();
            throw rpas::IllegalParse("Necessary outputs are missing");
         }


		  rpas::Intersection outputInt = _lhs[_args[RdfSvmExprArgs::SVM_OUTPUT]]->intersection();
		 
		 // check the input measures intersection
		 rpas::Intersection runMaskInt = _rhs[_args[RdfSvmExprArgs::RUN_MASK]]->intersection();
         for (int i=RdfSvmExprArgs::SVM_TYPE;i<RdfSvmExprArgs::DATA_HIER;i++)
		 {
			 rpas::Intersection tempInt = _rhs[_args[i]]->intersection();
			 if (runMaskInt !=tempInt)
			 {
               usage();
               throw rpas::IllegalParse("the RUM_MASK intersection should be same SVM Paramters measure");
			 }
		 }



		 rpas::Intersection targetInt = _rhs[_args[RdfSvmExprArgs::SVM_TARGET]]->intersection();
         for (int i=RdfSvmExprArgs::DATA_HIER +1;i<_totalArgNum -_RdfSvmExprFixedInputNum  ;i++)
		 {
			 rpas::Intersection tempInt = _rhs[_args[i]]->intersection();
			 if (targetInt !=tempInt)
			 {
               usage();
               throw rpas::IllegalParse("the SVM_TARGET intersection should be same SVM_FEATUREx measure");
			 }
		 }


		 if (outputInt != targetInt)
		 {
			    usage();
               throw rpas::IllegalParse("the SVM_TARGET intersection should be same SVM_OUTPUT measure");
		 }


		 rpas::Intersection outputGammaInt = _lhs[_args[RdfSvmExprArgs::SVM_OPTGAMMA]]->intersection();
		 rpas::Intersection outputConstInt = _lhs[_args[RdfSvmExprArgs::SVM_OPTCONST]]->intersection();

		 if (outputGammaInt !=runMaskInt ||
			 outputConstInt !=runMaskInt )
		 {
			 	 usage();
               throw rpas::IllegalParse("the SVM_OPTGAMMA,SVM_OPTCONST and SVM_RUNMASK intersection should be same ");
		 }

		 //targetInt.removeDimensionInfo(targetInt.retrieveDimensionAlong());

		
      }//virtual void specialExpressionValidateHook()


 void parseOutputMeasures(const std::vector<rpas::MeasureInstance>& atInstances)
    {


        _OutputArrayVec.resize(_RdfSvmExprOutputArgsNum );
        for (int i = 0; i < _RdfSvmExprOutputArgsNum ; ++i)
        {
			if (_args[i] !=-1)
			{
				rpas::VarParseNode* srcVpn = dynamic_cast<rpas::VarParseNode*>(_lhs[_args[i]]);
				if (srcVpn)
				{
					const rpas::MeasureWrapper& thisMeasWrapper =
							_measureWrapperManager.findWrapper(srcVpn->buildExternalName());
					_OutputArrayVec[i] = new rpas::SimpleArrayWrapper(atInstances[_args[i]],
							rpas::SimpleArrayWrapper::WRITE, _ExpressionEvalMode);
				}
				else
				{
					throw rpas::UnrecoverableExpressionLibraryError
							("In PopDataCount Special Expression: Syntax Error");
				}
			}
        }
    }


 void parseInputMeasures(const std::vector<rpas::MeasureInstance>& atInstances)
    {
		_InputArrays.resize(_totalArgNum);
        for (int i = _RdfSvmExprOutputArgsNum ; i < _totalArgNum ; i++)
        {
            rpas::logger() << rpas::debug << "Parsing argument " << i << "." << rpas::endl;
            int pos = _args[i];

            if (pos != -1)
            {
                rpas::VarParseNode* srcVpn = dynamic_cast<rpas::VarParseNode*>(_rhs[pos]);
                if (srcVpn)
                {
                    const rpas::MeasureWrapper& thisMeasWrapper =
                            _measureWrapperManager.findWrapper(srcVpn->buildExternalName());
                    rpas::MeasureInstance thisMeasInstance =
                            thisMeasWrapper.measureInstance(thisMeasWrapper.intersection());
                    //_InputArrays[i] = thisMeasInstance.array();
					_InputArrays[i] =srcVpn->getDefaultEvalMeasureInstance().array();
                     _InputArrays[i].mode(rpas::Ali::Read);
                }
                else
                {
					 // hier name may be a scalar integer
					rpas::String hierName ;
					rpas::ConstParseNode* cpnHierName = dynamic_cast<rpas::ConstParseNode *>(_rhs[pos]);
					if (cpnHierName)
					{
					   
					   hierName = cpnHierName->getStringVal();
					}

					if (i == RdfSvmExprArgs::DATA_HIER )
					{
						_dataHier = hierName;
					}

                }
            }
            rpas::logger() << rpas::debug << "Finished parsing argument " << i << "." << rpas::endl;
        }
   }

  void buildArrayMap(const rpas::DimensionSpace& fromSpace,
                            const rpas::DimensionSpace& toSpace,
                            rpas::ArrayMap* &flagmap)
{
    rpas::MeasureStore& measStore = rpas::MeasureStore::current();
    flagmap = new rpas::ArrayMap(fromSpace, toSpace, measStore);

}


 void initialization(const std::vector<rpas::MeasureInstance>& atInstances)
   {


	           //LHS 
        parseOutputMeasures(atInstances);

        //Open RHS Args in read mode
        parseInputMeasures(atInstances);


        buildArrayMap(_InputArrays[RdfSvmExprArgs::RUN_MASK].dimSpace(),_InputArrays[RdfSvmExprArgs::SVM_TARGET].dimSpace(),_runMaskToTarget);

			
		rpas::Intersection targetIntx=_rhs[_args[RdfSvmExprArgs::SVM_TARGET]]->intersection();
		rpas::Intersection runMaskIntx=_rhs[_args[RdfSvmExprArgs::RUN_MASK]]->intersection();
		targetIntx.removeDimensionInfo(targetIntx.retrieveDimensionAlong(_dataHier));
		if (targetIntx !=runMaskIntx)
		{

			usage();
			throw rpas::IllegalParse("RUN_MASK's intersection  should be SVM_TARGET intersection with removal of DATA_HIER ");
		}

		rpas::DimensionSpace dsMap =_InputArrays[RdfSvmExprArgs::SVM_TARGET].dimSpace();
		rpas::String dataDimName;
		_dataDimIndexInTarget = GetDimInfo(dsMap, _dataHier,_dataDimSize ,dataDimName);

   }

 bool autoConformRHS() const { 
	 return false; 
 }

virtual void specialExpressionEvalHook(
       const std::vector<rpas::MeasureInstance>& atInstances)
   {
		  try{
	   
		   // Set up the arrays
		   initialization(atInstances);
		   
		   rpas::Array maskArray = _InputArrays[RdfSvmExprArgs::RUN_MASK];

		   rpas::ArrayKey runMaskKey =  _InputArrays[RdfSvmExprArgs::RUN_MASK].emptyKey();
	   
		   rpas::ConstLogicalIterator runMaskIt = maskArray.constLogicalBegin();

		   rpas::ArrayCell ac;
		  for (; runMaskIt != _InputArrays[RdfSvmExprArgs::RUN_MASK].constLogicalEnd(); ++runMaskIt)
		  {
			 runMaskKey = runMaskIt.key();
			 maskArray.get(runMaskKey,ac);
			 if (ac.toBoolean()==false)
			 {
				 continue;
			 }
			
			  _InputArrays[RdfSvmExprArgs::SVM_TYPE].get(runMaskKey,ac);
			  int svm_type=ac.toInt();
			  if (svm_type != EPSILON_SVR &&svm_type !=NU_SVR)
			  {
				  rpas::logger() << rpas::debug << "wrong svm_type" << rpas::endl;
				  continue;
			  }

			  _InputArrays[RdfSvmExprArgs::SVM_KERNEL_TYPE].get(runMaskKey,ac);
			  int svm_kenel_type=ac.toInt();
			  if (svm_kenel_type != LINEAR &&svm_kenel_type !=POLY &&
				  svm_kenel_type!=RBF&&svm_kenel_type!=SIGMOID&&!svm_kenel_type !=PRECOMPUTED)
			  {
				  rpas::logger() << rpas::debug << "wrong svm_keenl_type" << rpas::endl;
				  continue;
			  }

			  _InputArrays[RdfSvmExprArgs::SVM_DEGREE].get(runMaskKey,ac);
			  int svm_degree=ac.toInt();
			  if (svm_degree<0)
			  {
				  rpas::logger() << rpas::debug << "wrong svm degree value" << rpas::endl;
				  continue;
			  }

			  _InputArrays[RdfSvmExprArgs::SVM_GAMMA].get(runMaskKey,ac);
			  double svm_gamma=ac.toNumeric();
			  if (svm_gamma<0)
			  {
				  rpas::logger() << rpas::debug << "wrong svm gamma value" << rpas::endl;
				  continue;
			  }

			  _InputArrays[RdfSvmExprArgs::SVM_COEF].get(runMaskKey,ac);
			  double svm_coeff=ac.toNumeric();

			  _InputArrays[RdfSvmExprArgs::SVM_CACHESIZE].get(runMaskKey,ac);
			  int svm_cache_size=ac.toInt();
			  if (svm_cache_size<0)
			  {
				  rpas::logger() << rpas::debug << "wrong svm cache size value" << rpas::endl;
				  continue;
			  }

			  _InputArrays[RdfSvmExprArgs::SVM_EPS].get(runMaskKey,ac);
			  double svm_eps=ac.toNumeric();
			  if (svm_eps<0)
			  {
				  rpas::logger() << rpas::debug << "wrong svm eps value" << rpas::endl;
				  continue;
			  }

			  _InputArrays[RdfSvmExprArgs::SVM_CONST].get(runMaskKey,ac);
			  double svm_const=ac.toNumeric();
			  if (svm_const<=0 &&((svm_type == C_SVC ||
								svm_type == EPSILON_SVR ||
								svm_type == NU_SVR)))
			  {
				  rpas::logger() << rpas::debug << "wrong svm constant value" << rpas::endl;
				  continue;
			  }

			  _InputArrays[RdfSvmExprArgs::SVM_NU].get(runMaskKey,ac);
			  double svm_nu=ac.toNumeric();
			  if ((svm_nu<=0|| svm_nu > 1) &&((svm_type == NU_SVC ||
								svm_type == ONE_CLASS ||
								svm_type == NU_SVR)))
			  {
				  rpas::logger() << rpas::debug << "wrong svm nu value" << rpas::endl;
				  continue;
			  }

			  _InputArrays[RdfSvmExprArgs::SVM_EPSION].get(runMaskKey,ac);
			  double svm_epsilon=ac.toNumeric();
			  if ((svm_epsilon<0) &&((svm_type == EPSILON_SVR)))
			  {
				  rpas::logger() << rpas::debug << "wrong svm epsilon value" << rpas::endl;
				  continue;
			  }

			  _InputArrays[RdfSvmExprArgs::SVM_SHRINKING].get(runMaskKey,ac);
			  int svm_shrink=ac.toInt();
			  if (svm_shrink!=0&&svm_shrink!=1)
			  {
				  rpas::logger() << rpas::debug << "wrong svm shrink value" << rpas::endl;
				  continue;
			  }

			  _InputArrays[RdfSvmExprArgs::SVM_ALPHA].get(runMaskKey,ac);
			  double svm_alpha=ac.toNumeric();
			  if (svm_alpha<0)
			  {
				  rpas::logger() << rpas::debug << "wrong svm alpha value" << rpas::endl;
				  continue;
			  }

			  _InputArrays[RdfSvmExprArgs::SVM_TRAIN_BEGIN].get(runMaskKey,ac);
			  int svm_train_begin_index=ac.toInt();

			  _InputArrays[RdfSvmExprArgs::SVM_TRAIN_END].get(runMaskKey,ac);
			  int svm_train_end_index=ac.toInt();

			  if(svm_train_end_index<svm_train_begin_index)
			  {
				 rpas::logger() << rpas::debug << "SVM_TRAIN_BEGIN should be less than SVM_TRAIN_END " << rpas::endl;
				  continue;
			  }

			  _InputArrays[RdfSvmExprArgs::SVM_PREDICT_BEGIN].get(runMaskKey,ac);
			  int svm_predict_begin_index=ac.toInt();

			  _InputArrays[RdfSvmExprArgs::SVM_PREDICT_END].get(runMaskKey,ac);
			  int svm_predict_end_index=ac.toInt();

			  if(svm_predict_end_index<svm_predict_begin_index)
			  {
				 rpas::logger() << rpas::debug << "SVM_PREDICT_BEGIN should be less than SVM_PREDICT_END " << rpas::endl;
				  continue;
			  }

			  _forecastLength=svm_predict_end_index - svm_predict_begin_index +1; 

			  _InputArrays[RdfSvmExprArgs::SVM_OPERATION].get(runMaskKey,ac);
			  int svm_operation=ac.toInt();		
			  if (svm_operation !=PREDICT && svm_operation !=OPT_PARAMETER
				                          && svm_operation !=PRE_OPT)
			  {
				  rpas::logger() << rpas::debug << "SVM_OPERATION should be one of {PREDICT,OPT_PARAMETER,BOTH} " << rpas::endl;
				  continue;
			  }

			  _InputArrays[RdfSvmExprArgs::SVM_NFOLDER].get(runMaskKey,ac);
			  int svm_nfolder=ac.toInt();
			  if (svm_nfolder<1&&(svm_operation ==OPT_PARAMETER ||svm_operation==PRE_OPT ))
			  {
				  rpas::logger() << rpas::debug << "the SVM_NFOLDER value should be greater than 0 for optimization mode" << rpas::endl;
				  continue;
			  }


			  struct svm_parameter param;	
			  param.degree=svm_degree;
			  param.svm_type=svm_type;
			  param.kernel_type=svm_kenel_type;
			  param.gamma = svm_gamma;
			  param.coef0=svm_coeff;
			  param.nu=svm_nu;
			  param.cache_size=svm_cache_size;
			  param.C = svm_const;
			  param.eps = svm_eps;
			  param.p=svm_epsilon;
			  param.shrinking = svm_shrink;
			  param.probability = 0;


			  
			   

			  rpas::Array targetArray =  _InputArrays[RdfSvmExprArgs::SVM_TARGET];
			  rpas::ArrayKey targetKey = targetArray.emptyKey();
			  int endIndex=0;
			  if (svm_train_end_index < _dataDimSize)
			  {
				  endIndex = svm_train_end_index;
			  }
			  else
			  {
				  endIndex = _dataDimSize;
			  }


			  bool hasmap = _runMaskToTarget->first(runMaskKey, targetKey);

			  double targetNaValue = _InputArrays[RdfSvmExprArgs::SVM_TARGET].naValue();

			  int trainLength=svm_train_end_index - svm_train_begin_index +1;
			  std::vector<double> trainDataVec;
			  std::vector<double> predictDataVec;
			  readData(targetKey,_dataDimIndexInTarget,svm_train_begin_index,svm_train_end_index,_InputArrays[RdfSvmExprArgs::SVM_TARGET],trainDataVec);

			  bool hasVaildData=false;
			  for (int i=0;i<trainDataVec.size();i++)
			  {
				  if (trainDataVec[i] != targetNaValue)
				  {
					  hasVaildData=true;
					  break;
				  }
			  }


			  if (hasVaildData==false)
			  {
				  continue;
			  }

			  std::vector<std::vector<double> > featuresTrainVec;
			  featuresTrainVec.resize(_dynamicFeaturesNum);

			  std::vector<std::vector<double> > featuresPredictVec;
			  featuresPredictVec.resize(_dynamicFeaturesNum);

			  _feaureNaValue=_InputArrays[RdfSvmExprArgs::DATA_HIER+1].naValue();


			  
			  for(int i =1;i<=_dynamicFeaturesNum;i++)
			  {
				  std::vector<double> trainFeatureVec;
				  readData(targetKey,_dataDimIndexInTarget,svm_train_begin_index,svm_train_end_index,_InputArrays[RdfSvmExprArgs::DATA_HIER+i],trainFeatureVec);
				  featuresTrainVec[i-1]=trainFeatureVec;
				  std::vector<double> predictFeatureVec;
				  readData(targetKey,_dataDimIndexInTarget,svm_predict_begin_index,svm_predict_end_index,_InputArrays[RdfSvmExprArgs::DATA_HIER+i],predictFeatureVec);
				  featuresPredictVec[i-1]=predictFeatureVec;
				  
			  }

			 /*
			     assume all the feature or target value has been scaled into [-1.0,10] or [0,1] before call into this special expression. 
				 the scaling could be easiliy done via Rpas Rules.
			 */
			  
			  if (svm_operation==OPT_PARAMETER ||
				  svm_operation==PRE_OPT)
			  {
				  optmizeParameter(param,trainDataVec, featuresTrainVec,_forecastLength,svm_nfolder,svm_alpha);
			  }

			  if (svm_operation==OPT_PARAMETER ||
				  svm_operation==PRE_OPT)
			  {
				  rpas::Array outpuGammatArray = _OutputArrayVec[RdfSvmExprArgs::SVM_OPTGAMMA]->getDataArray();
				  rpas::Array outpuConstArray = _OutputArrayVec[RdfSvmExprArgs::SVM_OPTCONST]->getDataArray();
				  outpuGammatArray.set(runMaskKey,param.gamma);
				  outpuConstArray.set(runMaskKey,param.C);
				  
			  }

			  if (svm_operation==OPT_PARAMETER )
			  {
				svm_free_and_destroy_model(&_model);
				featuresTrainVec.clear();
				predictDataVec.clear();
				continue;
			  }


			  buildProblem(_prob,trainDataVec,featuresTrainVec);
			  
			  _model=Malloc(svm_model,1);
			  _model = svm_train(&_prob,&param);
			  //svm_save_model("savedModel.txt",_model);
			  rpas::Array outputArray = _OutputArrayVec[RdfSvmExprArgs::SVM_OUTPUT]->getDataArray();
			  rpas::ArrayKey outputKey =  outputArray.emptyKey();

			  std::vector<double> outputVec;
			  predictSVMValue(_model,svm_predict_end_index -svm_predict_begin_index +1,featuresPredictVec,outputVec);
			  

				//outputKey.set(_dataDimIndexInTarget,svm_predict_begin_index+i);
				//outputArray.set(outputKey,predictValue);
			  writeOutput(runMaskKey,outputArray,svm_predict_begin_index,outputVec);

			  svm_free_and_destroy_model(&_model);
			  featuresTrainVec.clear();
			  predictDataVec.clear();
	          free(_prob.y);
			  for (int i=0;i<_prob.l;i++)
			  {
				  if (_prob.x[i] !=0)
				  {
					  free(_prob.x[i]);
				  }
			  }
			  free(_prob.x);


		  }
	   }
	   catch(...)
	   {
		   cleanup();
		   throw;
	   }
   }

   void optmizeParameter(struct svm_parameter& para,const std::vector<double> & trainDataVec, const std::vector<std::vector<double> >& trainFeatureVec,int folderWidth,int backcastSteps,double alpha)
   {
	   std::vector<double> CList, GList;
		double baseNum = 2.0;
		for(double j = -5; j <= 15; j += 2) //-5 and 15
		CList.push_back(pow(baseNum,j));
		for(double j = -15; j <= 3; j += 2) //-15 and 3
		 GList.push_back(pow(baseNum,j));

		if (alpha<=0)
		{
			folderWidth /trainDataVec.size();
		}
		double minMse=std::numeric_limits<double>::max();
		double bestC=0;
		double bestGamma=0;

		for(int i= 0; i<CList.size();i++) //for all C's
		{
			double C = CList[i];
			for(int j=0;j<GList.size();j++) //for all gamma's
			{
				double gamma = GList[j];
				para.C = C;
				para.gamma = gamma;

				double sumError=0;

				for (int k=0;k<backcastSteps;k++)
				{

					if (trainDataVec.size() - (k+1)*folderWidth < folderWidth)
					{
						break;
					}

					svm_problem currentProblem;
					std::vector<double> curTrainDataVec(trainDataVec.size() -(k+1)*folderWidth);
					std::vector<double> curPredictDataVec(folderWidth);
					std::copy(trainDataVec.begin(),trainDataVec.end() -(k+1)*folderWidth,curTrainDataVec.begin());
					std::copy(trainDataVec.end() -(k+1)*folderWidth,trainDataVec.end() -k*folderWidth,curPredictDataVec.begin());

					std::vector<std::vector<double> > CurTrainFeatureVec(trainFeatureVec.size());
					std::vector<std::vector<double> > CurPredictFeatureVec(trainFeatureVec.size());

					for (int f=0;f<trainFeatureVec.size();f++)
					{
						std::vector<double>tempTrainFeatureVec(trainDataVec.size() -(k+1)*folderWidth);
						std::vector<double>tempPredictFeatureVec(folderWidth);
						std::copy(trainFeatureVec[f].begin(),trainFeatureVec[f].end() -(k+1)*folderWidth,tempTrainFeatureVec.begin());
						std::copy(trainFeatureVec[f].end() -(k+1)*folderWidth,trainFeatureVec[f].end() -k*folderWidth,tempPredictFeatureVec.begin());
						CurTrainFeatureVec[f]=tempTrainFeatureVec;
						CurPredictFeatureVec[f]=tempPredictFeatureVec;
					}

			
					
					buildProblem(currentProblem,curTrainDataVec,CurTrainFeatureVec);
			  
					svm_model * curModel=Malloc(svm_model,1);
					curModel = svm_train(&currentProblem,&para);
					std::vector<double> outputVec;
					//outputVec.resize(folderWidth);
					predictSVMValue(curModel,folderWidth,CurPredictFeatureVec,outputVec);

					double currMse,CurrR2;
					svmMseAndR2(curPredictDataVec,outputVec,currMse,CurrR2);
					sumError +=pow(alpha,i)*currMse;


					svm_free_and_destroy_model(&curModel);
			        free(currentProblem.y);
					for (int i=0;i<currentProblem.l;i++)
					{
						if (currentProblem.x[i] !=0)
						{
							free(currentProblem.x[i]);
						}
					}
					free(currentProblem.x);

				}

				if (sumError <minMse)
				{
					bestC=C;
					bestGamma =gamma;
				}

			}
		}


		para.C = bestC;
		para.gamma = bestGamma;
   }


   void writeOutput(rpas::ArrayKey ak,rpas::Array outputArray,int startIndex,std::vector<double>outputValuesVec)
   {
	   rpas::ArrayKey outputKey = outputArray.emptyKey();
	   bool hasmap = _runMaskToTarget->first(ak, outputKey);
	   for(int i=0;i<outputValuesVec.size();i++)
	   {
		   	outputKey.set(_dataDimIndexInTarget,startIndex+i);
			double predictValue=outputValuesVec[i];
			outputArray.set(outputKey,predictValue);

	   }
	   
   }

   void svmMseAndR2(const std::vector<double> & y,const std::vector<double> & x,double& mse,double& r2)
   {
	   if (y.size() !=x.size())
	   {
		    rpas::logger() << rpas::debug << "the y.size should be equal to predict.size " << rpas::endl;
			return;
	   }

	   double x_y=0,xy=0,sum_x=0,sum_y=0,sum_xx=0,xx=0,yy=0,sum_yy=0;
	   long l=x.size();
	   for (int i=0;i<x.size();i++)
	   {
		   x_y +=(x[i]-y[i])*(x[i]-y[i]);
		   xy +=x[i]*y[i];
		   sum_x +=x[i];
		   sum_y +=y[i];
		   sum_xx +=x[i]*x[i];
		   sum_yy +=y[i]*y[i];
	   }
	   mse = x_y/l;
	   r2=(xy/l-sum_x*sum_y)*(xy/l-sum_x*sum_y)/((sum_xx/l -sum_x*sum_x)*(sum_yy/l - sum_y*sum_y));
   }

	void predictSVMValue(svm_model * model,int forecastLength,const std::vector<std::vector<double> >& featuresPredictVec,std::vector<double>& outputVec)
	 {
	  		  for(int i=0;i<forecastLength;i++)
			  {
				svm_node * x = Malloc(struct svm_node,_dynamicFeaturesNum+1 );
				int validFeatureCount=0;
				for(int fIndex=0;fIndex<featuresPredictVec.size();fIndex++)
				{
					if(featuresPredictVec[fIndex][i] != _feaureNaValue)
					{
						x[validFeatureCount].index=fIndex;
						x[validFeatureCount].value = featuresPredictVec[fIndex][i];
						validFeatureCount++;
					}

				}

				x[validFeatureCount].index = -1;

				double predictValue=svm_predict(model,x);

				//outputKey.set(_dataDimIndexInTarget,svm_predict_begin_index+i);
				//outputArray.set(outputKey,predictValue);
				outputVec.insert(outputVec.end(),predictValue);
				free(x);
			  }
  }

   void readData(rpas::ArrayKey ak,int dimIndex,int fromIndex,int toIndex,rpas::Array ar,std::vector<double> & toVec)
   {
	   rpas::ArrayKey currentKey= ak;

	   rpas::ArrayCell ac;
	   for (int i=fromIndex;i<=toIndex;i++)
	   {
		   currentKey.set(dimIndex,i);
		   ar.get(currentKey,ac);
		   toVec.insert(toVec.end(),ac.toNumeric());
	   }
   }

   void buildProblem(struct svm_problem& sp,const std::vector<double>& data,const std::vector<std::vector<double> >& features)
   {
	   int count=0;
	   for(int i=0;i<data.size();i++)
	   {
		   if (data[i] !=0)
		   {
			   count++;
		   }
	   }

	   sp.l=count;
	   
	   sp.y = Malloc(double,sp.l);
	   sp.x = Malloc(struct svm_node *,sp.l);
	   
	   int currentYIndex=0;

	   for(int i=0;i<data.size();i++)
	   {
		   if (data[i]==0)
		   {
				continue;
		   }
		   else
		   {
				sp.y[currentYIndex]=data[i];
				currentYIndex++;
		   }
		   //count how many valid feature for current y
		   int validFeatureCount=0;
		   for (int f=0;f<_dynamicFeaturesNum;f++)
		   {
			   if(features[f][i]!=_feaureNaValue)
			   {
				   validFeatureCount++;
			   }
		   }

		   svm_node * x_space = Malloc(struct svm_node,validFeatureCount+1 );
		   sp.x[i] = &x_space[0];
		   // put the current feature into the x
		   validFeatureCount=0;
		   for (int f=0;f<_dynamicFeaturesNum;f++)
		   {
			   if(features[f][i]!=_feaureNaValue)
			   {
				   x_space[validFeatureCount].index=f;
				   x_space[validFeatureCount].value=features[f][i];
				   validFeatureCount++;
			   }
		   }

		   x_space[validFeatureCount].index=-1;

	   }
		 



   }

       // Overload for expression incremental evaluation
   virtual void specialExpressionEvalIncrementalHook(
       const std::vector<rpas::MeasureInstance>& atInstances)
   {
      initialization(atInstances);
	  _ExpressionEvalMode = rpas::ExpressionEvalMethod::incremental();
      specialExpressionEvalHook(atInstances);
   }

 
	  //-----------general function--------------------------------

      int GetDimInfo(const rpas::Intersection& thisInt,
					 const rpas::String& hierName,
					 int& dimSize,
					 rpas::String& dimName)
      {
		const rpas::DimensionSpace& dimSpace = thisInt.getDimensionSpace();
		int dimnumber = dimSpace.dimensionality(); 

		rpas::MeasureStore& ms = rpas::MeasureStore::current();
		for (int i = 0; i < dimnumber-1; ++i)
		{
			rpas::Dimension* dim = dimSpace.dimension(i);
			rpas::DimensionInfo& diminfo = ms.getDimensionInfo(dim->name());
			if (hierName.equalsIgnoreCase(diminfo.hierarchy()))
			{
				dimName = dim->name();
				dimSize = dim->size();
				return i;
			}
		}
		return -1;
      }

	  int GetDimInfo(const rpas::DimensionSpace& dimSpace,
                           const rpas::String& hierName,
                           int& dimSize,
                           rpas::String& dimName)	
	  {
        int dimnumber = dimSpace.dimensionality(); 

        rpas::MeasureStore& ms = rpas::MeasureStore::current();
        for (int i = 0; i < dimnumber-1; ++i)
        {
            rpas::Dimension* dim = dimSpace.dimension(i);
            rpas::DimensionInfo& diminfo = ms.getDimensionInfo(dim->name());
            if (hierName.equalsIgnoreCase(diminfo.hierarchy()))
            {
                dimName = dim->name();
                dimSize = dim->size();
                return i;
            }
        }
        return -1;
	}

      void usage()
      {
         rpas::logger() << rpas::error << rpas::String(" Usage of RdfSvmExpr")<<rpas::endl;
         int i;
         bool firstInput = false;
         bool prevOut = true;

         for (i =0; i<_totalArgNum  ;i++)
         {
            apps::ExprParaInfo para;
            _args.paraInfoLookUp(i,para);

            if (i>0 && i<_totalArgNum -1)
               if (!(prevOut && para._input))
                  rpas::logger()<<rpas::error<<","<<rpas::endl;
               else
                  rpas::logger()<<rpas::error<<rpas::endl;
            else
               rpas::logger()<<rpas::error<<rpas::endl;


            if (!para._input && !firstInput)
            {
               rpas::logger()<<rpas::error<<"<-RdfSvmExpr("<<rpas::endl;
               firstInput = true;
            }

         }
         rpas::logger()<<rpas::error<<")"<<rpas::endl;
      }

      //-----------general function--------------------------------


	};//class RdfSvmExpr
	

	static int _RdfSvmExprId = 
		rpas::Parser::registerSpecialExpressionFactory(
		new rpas::SimpleSpecialExpressionFactory<RdfSvmExpr>(
		"RdfSvmExpr"));


}//end namespace rdf

