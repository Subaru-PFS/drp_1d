#include <RedshiftLibrary/common/quicksort.h>
// function template
template <typename T>
CFindExtrema<T>::CFindExtrema()
{
}
template <typename T>
CFindExtrema<T>::~CFindExtrema()
{
}

//A great part of below code is meant to disappear after merging fix_issue_5619, since the call to CExtreme::Find will include the refine-search 
template <class T>
void CFindExtrema<T>::callFind(T chisquareResult ){
            Int32 extremumCount = 10;
            if (chisquareResult->Redshifts.size() > extremumCount)
            {
                TPointList extremumList;
                TFloat64Range redshiftsRange(
                 chisquareResult->Redshifts[0],
                 chisquareResult->Redshifts[chisquareResult->Redshifts.size() - 1]);
                CExtremum extremum(redshiftsRange, extremumCount, true);
                extremum.Find(chisquareResult->Redshifts, chisquareResult->ChiSquare, extremumList);

            //*
            // Refine Extremum with a second maximum search around the z candidates:
            // This corresponds to the finer xcorrelation in EZ Pandora (in
            // standard_DP fctn in SolveKernel.py)
                Float64 radius = 0.001;
                for (Int32 i = 0; i < extremumList.size(); i++)
                {
                    Float64 x = extremumList[i].X;
                    Float64 left_border = max(redshiftsRange.GetBegin(), x - radius);
                    Float64 right_border = min(redshiftsRange.GetEnd(), x + radius);

                    TPointList extremumListFine;
                    TFloat64Range rangeFine = TFloat64Range(left_border, right_border);
                    CExtremum extremumFine(rangeFine, 1, true);
                    extremumFine.Find(chisquareResult->Redshifts, chisquareResult->ChiSquare,
                              extremumListFine);
                    if (extremumListFine.size() > 0)
                    {
                    extremumList[i] = extremumListFine[0];
                    }
                }
                // store extrema results
                chisquareResult->Extrema.resize(extremumCount);
                for (Int32 i = 0; i < extremumList.size(); i++)
                {

                    chisquareResult->Extrema[i] = extremumList[i].X;
                }

            } else
            {
                // store extrema results
                chisquareResult->Extrema.resize(chisquareResult->Redshifts.size());
                TFloat64List tmpX;
                TFloat64List tmpY;
                for (Int32 i = 0; i < chisquareResult->Redshifts.size(); i++)
                {
                    tmpX.push_back(chisquareResult->Redshifts[i]);
                    tmpY.push_back(chisquareResult->ChiSquare[i]);
                }
                // sort the results by merit
                NSEpic::CQuickSort<Float64> sort;
                vector<Int32> sortedIndexes(chisquareResult->Redshifts.size());
                sort.SortIndexes(tmpY.data(), sortedIndexes.data(),
                         sortedIndexes.size());
                for (Int32 i = 0; i < chisquareResult->Redshifts.size(); i++)
                {
                    chisquareResult->Extrema[i] = tmpX[sortedIndexes[i]];
                }
            }
        return;
}

