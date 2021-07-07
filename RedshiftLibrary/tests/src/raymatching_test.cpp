#include "RedshiftLibrary/operator/raymatching.h"
#include "RedshiftLibrary/ray/catalog.h"
#include "RedshiftLibrary/ray/ray.h"


#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;
BOOST_AUTO_TEST_SUITE(test_raymatching)


Float64 precision = 1e-12;
BOOST_AUTO_TEST_CASE(Compute){

  CRayCatalog detectedRayCatalog;
  std::shared_ptr<CLineProfile> profilesym{std::make_shared<CLineProfileSYM>()};
  detectedRayCatalog.Add(CRay ( CRay("Ray1",20, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  detectedRayCatalog.Add(CRay ( CRay("Ray2",40, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));

  detectedRayCatalog.Add(CRay ( CRay("Ray3",250, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  detectedRayCatalog.Add(CRay ( CRay("Ray4",800, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));


  CRayCatalog restRayCatalog;
  restRayCatalog.Add(CRay ( CRay("Ray2_1",10, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  restRayCatalog.Add(CRay ( CRay("Ray2_2",20, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  restRayCatalog.Add(CRay ( CRay("Ray2_3",125, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  restRayCatalog.Add(CRay ( CRay("Ray2_4",400, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));

  CRayMatching lineMathing = CRayMatching();
  TFloat64Range redshiftRange = TFloat64Range(0,5);
  Int32 nThreshold =4;
  Float64 tol =0.05;
  Int32 typeFilter = CRay::nType_Emission;
  Int32 detectedForceFilter = 2;
  Int32 restRorceFilter = 2;

  std::shared_ptr<CRayMatchingResult> res = lineMathing.Compute(detectedRayCatalog, restRayCatalog, redshiftRange, nThreshold, tol, typeFilter, detectedForceFilter, restRorceFilter);
  BOOST_CHECK_EQUAL(res->SolutionSetList.size(), 1);

  //BOOST_CHECK_EQUAL(res->SolutionSetList[0].size(), 1);
  BOOST_CHECK_CLOSE(res->SolutionSetList[0][0].Redshift,1, precision);


  detectedRayCatalog = CRayCatalog();
  detectedRayCatalog.Add(CRay ( CRay("Ray1",800, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));


  restRayCatalog= CRayCatalog();
  restRayCatalog.Add(CRay ( CRay("Ray2_1",10, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  restRayCatalog.Add(CRay ( CRay("Ray2_2",20, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  restRayCatalog.Add(CRay ( CRay("Ray2_3",125, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  restRayCatalog.Add(CRay ( CRay("Ray2_4",400, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  nThreshold =1;
  res = lineMathing.Compute(detectedRayCatalog, restRayCatalog, redshiftRange, nThreshold, tol, typeFilter, detectedForceFilter, restRorceFilter);
  BOOST_CHECK_EQUAL(res->SolutionSetList.size(), 1);

  //BOOST_CHECK_EQUAL(res->SolutionSetList[0].size(), 1);
  BOOST_CHECK_CLOSE(res->SolutionSetList[0][0].Redshift,1, precision);


  detectedRayCatalog = CRayCatalog();
  detectedRayCatalog.Add(CRay ( CRay("Ray1",20, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  detectedRayCatalog.Add(CRay ( CRay("Ray2",20, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  detectedRayCatalog.Add(CRay ( CRay("Ray3",40, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  detectedRayCatalog.Add(CRay ( CRay("Ray4",80, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  detectedRayCatalog.Add(CRay ( CRay("Ray5",160, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  detectedRayCatalog.Add(CRay ( CRay("Ray6",320, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));


  restRayCatalog= CRayCatalog();
  restRayCatalog.Add(CRay ( CRay("Ray2_1",10, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  restRayCatalog.Add(CRay ( CRay("Ray2_2",10, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  restRayCatalog.Add(CRay ( CRay("Ray2_3",20, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  restRayCatalog.Add(CRay ( CRay("Ray2_4",40, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  restRayCatalog.Add(CRay ( CRay("Ray2_5",80, CRay::nType_Emission, profilesym, 2, 1, 4, 5.6)));
  nThreshold =3;
  redshiftRange = TFloat64Range(0,10);
  res = lineMathing.Compute(detectedRayCatalog, restRayCatalog, redshiftRange, nThreshold, tol, typeFilter, detectedForceFilter, restRorceFilter);
  BOOST_CHECK_EQUAL(res->SolutionSetList.size(), 3);

  //BOOST_CHECK_EQUAL(res->SolutionSetList[0].size(), 1);
  BOOST_CHECK_CLOSE(res->SolutionSetList[0][0].Redshift,1, precision);
  BOOST_CHECK_CLOSE(res->SolutionSetList[1][0].Redshift,3, precision);
  BOOST_CHECK_CLOSE(res->SolutionSetList[2][0].Redshift,7, precision);
}

BOOST_AUTO_TEST_SUITE_END()
