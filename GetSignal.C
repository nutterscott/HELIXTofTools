// methods to extract arithmetic and geometric means from paddle ends and full paddles
//in this file:
//  double getPadEndMeanSignal(HLX_Stg2_TOFPaddleEndHit endHit, bool useArith = true) const  //default to using arithmetic mean
//  double getPadMeanSignal(HLX_Stg2_TOFHit* padHit, bool usePadArith = false, bool useEndArith = true) const 
//  double getPadMeanArithSignal(HLX_Stg2_TOFHit* padHit, bool useEndArith = true) 
//  double getPadMeanGeomSignal(HLX_Stg2_TOFHit* padHit, bool useEndArith = true) 
//  double getPadEndArithMeanSignal(HLX_Stg2_TOFPaddleEndHit endHit)
//  double getPadEndGeomMeanSignal(HLX_Stg2_TOFPaddleEndHit endHit)

  
double getPadEndMeanSignal(HLX_Stg2_TOFPaddleEndHit endHit, bool useArith = true)  //default to using arithmetic mean
{
  double sumArithmetic{0};
  double sumGeometric{1};

  int    feeCount{0};
  for (auto &fee : endHit.fFEEHits) {
    if (fee.fStatusFlag == HLX_Stg2_TOFFEEHit::kBad) continue;
    if (fee.fStatusFlag == HLX_Stg2_TOFFEEHit::kCalibrationOnly) continue;
    feeCount++;
    sumArithmetic += fee.fSignalSum;
    sumGeometric  *= fee.fSignalSum;
  }
  if (useArith)
    return (feeCount==0)? 0 : sumArithmetic/feeCount;           //Arithmetic mean
  else
    return (feeCount==0)? 0 : pow(sumGeometric, 1.0/feeCount); //Geometric mean
}

double getPadEndArithMeanSignal(HLX_Stg2_TOFPaddleEndHit endHit) {return getPadEndMeanSignal(endHit, true);}
double getPadEndGeomMeanSignal (HLX_Stg2_TOFPaddleEndHit endHit) {return getPadEndMeanSignal(endHit, false);}

double getPadMeanSignal(HLX_Stg2_TOFHit* padHit, bool usePadArith = false, bool useEndArith = true) 
{
  HLX_Stg2_TOFPaddleEndHit tE     = padHit->getPaddleEndHit(HLX_Geometry_TOF::kEast);
  HLX_Stg2_TOFPaddleEndHit tW     = padHit->getPaddleEndHit(HLX_Geometry_TOF::kWest);
  
  double sE = (useEndArith)? getPadEndArithMeanSignal(tE) : getPadEndGeomMeanSignal(tE);
  double sW = (useEndArith)? getPadEndArithMeanSignal(tW) : getPadEndGeomMeanSignal(tW);
  if (usePadArith) return (sE+sW)/2.;
  else return (sE*sW>=0)? sqrt(sE*sW) : 0; 
}

double getPadMeanArithSignal(HLX_Stg2_TOFHit* padHit, bool useEndArith = true) {return getPadMeanSignal(padHit, true, useEndArith);}
double getPadMeanGeomSignal (HLX_Stg2_TOFHit* padHit, bool useEndArith = true) {return getPadMeanSignal(padHit, false, useEndArith);}

double getTBPadEndSignalSC(HLX_Stg2_TOFPaddleEndHit endHit)  //default to using arithmetic mean
{
  //arithmetic only for top and bottom
  double endSig = 0;
  int nFees     = endHit.getFEEHitCount();
  int nFEEgood  = 0;
  for (int j=0; j<nFees; j++)
  {
  	if ( endHit.getFEEHit(j).fStatusFlag == HLX_Stg2_TOFFEEHit::kGood ) {
  	  nFEEgood++;
  	  endSig += endHit.getFEEHit(j).fSignalSum;
  	}
  }
  if (nFEEgood == 0) return 0;

  endSig /= nFEEgood;  // now we have the ARITHMETIC average of the paddle end hit (at this point still messed up for the mpv[8] / mpv[0] misnormalization, this will be fixed
  
  return endSig;
}  

double getBoreEndSignalSC(HLX_Stg2_TOFPaddleEndHit endHit)  //default to using arithmetic mean
{
  //arithmetic only for top and bottom
  double endSig = 0;
  int nFees     = endHit.getFEEHitCount();

  int nFEEgood  = 0;
  for (int j=0; j<nFees; j++)
  {
  	if ( endHit.getFEEHit(j).fStatusFlag == HLX_Stg2_TOFFEEHit::kGood ) {
  	  nFEEgood++;
  	  endSig += endHit.getFEEHit(j).fSignalSum;
  	}
  }
  if (nFEEgood == 0) return 0;

  endSig /= nFEEgood;  // now we have the ARITHMETIC average of the paddle end hit (at this point still messed up for the mpv[8] / mpv[0] misnormalization, this will be fixed

  /*
  if (nFees != 4) return 0; 

	endSig = max(endHit.getFEEHit(0).fSignalSum+endHit.getFEEHit(1).fSignalSum,
     max(endHit.getFEEHit(1).fSignalSum+endHit.getFEEHit(2).fSignalSum,
	   endHit.getFEEHit(2).fSignalSum+endHit.getFEEHit(3).fSignalSum));
   endSig /= 2.;    // divide by 2 for average
  */
  return endSig;
}  
  
  
double getPadSignalSC(HLX_Stg2_TOFHit* padHit) 
{
  HLX_Stg2_TOFPaddleEndHit tE     = padHit->getPaddleEndHit(HLX_Geometry_TOF::kEast);
  HLX_Stg2_TOFPaddleEndHit tW     = padHit->getPaddleEndHit(HLX_Geometry_TOF::kWest);
  
  auto theGeom = HLX_Geometry::getInstance();
  double sE = 0, sW = 0;
  double padSig = 0;
  if(theGeom->fTOF.isBore(padHit->fTag.getPaddleID()))
  {  //bore paddle
    sE = getBoreEndSignalSC(tE);
    sW = getBoreEndSignalSC(tW);
    padSig = (sE+sW)/2.;
  } else { //TB paddle
    sE = getTBPadEndSignalSC(tE);
    sW = getTBPadEndSignalSC(tW);
    padSig = sqrt(sE*sW);
  }
  return padSig;
}
