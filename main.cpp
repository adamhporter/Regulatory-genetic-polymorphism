//
//  main.cpp
//  FitnessOverdominanceOptimalGenotypeByPopt
//
//  Created by adam porter on 4/4/19.
//  Copyright © 2019 adam porter. All rights reserved.
//
//	This program calculates the optimal 2-locus, 2-allele genotype
//	for a simple regulatory genetic pathway, when the
//	regulatory and structural sites of the transcription factor locus,
//	and the regulatory (cis) site of the expressed locus can vary.
//
//	Genotypes are represented by bitstrings of specifiable length.
//
//	This entails searching through all genotypic combinations, and for each:
//		-- determining the possible 2-locus, 2-allele genotypes in a population that contained that genotype
//		-- calculating the fitnesses of each of those genotypes, given an optimal fitness (Popt) and a fitness variance (omega)
//		-- maximizing the equation for population mean fitness to determine the optimal allele frequencies
//		-- sorting by population mean fitness to find the genotypic combination with the highest population mean fitness
//		-- saving that
//	This is repeated for for 0≤Popt≤1 in steps of 0.01, for a given omega.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <cstring>
#include <vector>
#include <climits>
#include <stdint.h>
#include <thread>
#include <utility>
#include <mutex>

std::recursive_mutex coutLock;//to make sure separate threads don't call std::cout simultaneously

class SimplestRegPathIndividual;
template <class T> inline int HammingDistance(T m, T n);
template <class T> inline T MIN(T x, T y);
template <class T> inline T MAX(T x, T y);

static int decimalDigitsToRound = 7;

static long double zero = (long double)0.0;
static long double half = (long double)0.5;
static long double one = (long double)1.0;
static long double two = (long double)2.0;
static long double four = (long double)4.0;
static long double ten = (long double)10.0;
//static uint64_t mask8bit[8]= {255,65280,16711680,4278190080,1095216660480,280375465082880,71776119061217280,18374686479671623680};
                        //mask8bit[]="0000000000000000000000000000000000000000000000000000000011111111";
                        //mask8bit[]="0000000000000000000000000000000000000000000000001111111100000000";
                        //mask8bit[]="0000000000000000000000000000000000000000111111110000000000000000";
                        //mask8bit[]="0000000000000000000000000000000011111111000000000000000000000000";
                        //mask8bit[]="0000000000000000000000001111111100000000000000000000000000000000";
                        //mask8bit[]="0000000000000000111111110000000000000000000000000000000000000000";
                        //mask8bit[]="0000000011111111000000000000000000000000000000000000000000000000";
                        //mask8bit[]="1111111100000000000000000000000000000000000000000000000000000000";
//static uint64_t mask10bit[6]= {1023,1047552,1072693248,1098437885952,1124800395214848,1151795604700004352};

template <class T> inline T ABS(T x)
	{return x<(T)0 ? (-x):x;}//ABSOLUTE VALUE TEMPLATE FUNCTION
template <class T> inline T MIN(T x, T y)
	{return x<y ? x:y;}
template <class T> inline T MAX(T x, T y)
	{return y<x ? x:y;}
template <class T> inline long MOD(T i, T j)
    {return long( i - j * long(i/j));}
template <class T> inline int SIGN(T x)
	{return x<0 ? (-1):1;}
template <class T> inline T SIGN(T a,T b)//Press et al. 1992 function
	{return (b >= zero ? ABS(a) : -ABS(a));}
template <class T> inline T BOUND(T x, T low, T high){
	x=MAX(x,low); return MIN(x,high);}
template <class T> inline void SWAP(T& x, T& y)
	{T temp = x; x = y; y = temp;}

template <class T> inline T ROUND(T n, const int decimalPlace){//n is the number: base 10
		long double val, nn = (long double)n, p=(long double)1.0,nsign=(long double)1.0;
		int i;
		if(nn<(long double)0.0){
			nsign=(long double)-1.0;
			nn *= nsign;}
		if(decimalPlace>=0){
				for(i=0;i<decimalPlace;++i){
					p *= (long double)10.0;}}
			else{
				for(i=0;i<-decimalPlace;++i){
					p /= (long double)10.0;}}
		nn *= p;
		val = round(nn)/p;
		if(val!=(long double)0.0){
			val*=nsign;}
		return (T)val;
		}//ROUND(T, const int)


void MonthTextToNum(std::string& textMo, int& numMo);// January || Jan -> 1
void MonthNumToText(std::string& textMo, int& numMo);// 1 -> January
void MonthNumToTextAbbreviated(std::string& textMo, int& numMo);// 1 -> Jan

class aTime{
    public:
        time_t firstTime_;//when class created
        time_t mostRecentTime_;//last call
        time_t elapsedSinceLast_;//when last call is replaced
        
    public:
        aTime(void){
            time(&firstTime_);//sets to current time
            mostRecentTime_=firstTime_;
            elapsedSinceLast_=0;
            }//end constructor
        
        aTime(aTime& t):firstTime_(t.firstTime_),
                mostRecentTime_(t.mostRecentTime_),
                elapsedSinceLast_(t.elapsedSinceLast_){}//copies

virtual ~aTime(void){}

virtual aTime* Duplicate(void){
            return new aTime(*this);}//virtual copies
        
    public:
        time_t first(void){return firstTime_;}
        time_t recent(void){return mostRecentTime_;}
        time_t elapsedSinceLast(void){return elapsedSinceLast_;}
        long double elapsed(void){//mostRecentTime_ to now, in seconds
            time_t now;
            time(&now);
            return difftime(now,mostRecentTime_);}
        
        void mark(void){//sets new most recent time
            time_t temp = mostRecentTime_;
            time(&mostRecentTime_);
            elapsedSinceLast_ = mostRecentTime_ - temp;
            }


        std::string* HMS(long double timediff, char sepchar=':'){
            //returns a string with hours:mins:secs only
            //used mostly to display differences in time
            long hours = (long) timediff/3600;
            long double td = timediff - (hours*(long double)3600);
            int mins = (int)td/60;
            td = td-((long double)mins*(long double)60);
            std::stringstream ssHr,ssMin,ssSec;
            ssHr<<hours<<sepchar;
            if(mins<10){ssMin<<'0';}
            ssMin<<mins;
            ssMin<<sepchar;
            int secs = (int) td;
            td = secs-td;//leaves decimal
            if(secs<10){ssSec<<'0';}
            ssSec<<secs;
            ssHr<<ssMin.str()<<ssSec.str();
            if(td >zero){
                std::stringstream ssTd;
                ssTd<<td;
                std::string tdStr=ssTd.str();
                while(tdStr[0]=='0'){
                    tdStr.erase(0,1);}
                ssHr<<tdStr;
                }
            std::string* hms = new std::string(ssHr.str());
            return hms;
            }//HMS

    std::string* HMS_elapsed(char sepchar=':'){
            return HMS(elapsed(),sepchar);}
            //returns a string with hours:mins:secs only
            //used for time since most recent marked time
        
        std::string* HMS_totalElapsed(char sepchar=':'){
            return HMS(totalElapsed(),sepchar);}
            //returns a string with hours:mins:secs only
            //used for time since first time
            
        time_t now(void){
            time_t moment; return time(&moment);}
            
            
        std::string* localStartTime(int *orderToReport=NULL){
            return localTimeAtMoment(firstTime_,orderToReport);}
            
            
        std::string* localTime(int *orderToReport=NULL){
            time_t moment = now();
            return localTimeAtMoment(moment,orderToReport);}
        
        
        std::string* localTimeAtMoment(time_t momentInTime, int *orderToReport=NULL){
                        //orderToReport elements: [0]=yr, [1]=mo, [2]=day, [3]=H, [4]=min, [5]=sec, [6]=use am/pm,
                        //[7]=include day of wk in text, [8]=use mo in text
                        //e.g. output:
                            //012345000 ==> 2013 10/05 15:31:45
                            //012345100 --> 2013 10/05 3:31:45 pm
                            //201345111 --> Sat October 5, 2013 3:31:45 pm
                        //order==null: default Sat Oct 5
            struct tm *moment = NULL;
            moment = localtime(&momentInTime);
            std::string timeStr = asctime(moment);
            std::string defaultDisplayTime(timeStr);//Wed Feb 13 17:17:11 2013
            while(defaultDisplayTime.back()=='\n'){
                defaultDisplayTime.pop_back();}
            if(orderToReport==NULL) {
                return new std::string(defaultDisplayTime);}
            std::string yr,d,h,m,s,day,month,colon(":"),displayTimeCopy(defaultDisplayTime);
            size_t spacepos=displayTimeCopy.find_first_of(" ");
            day=displayTimeCopy.assign(0,spacepos-1);
            displayTimeCopy.erase(0,spacepos);
            spacepos=displayTimeCopy.find_first_of(" ");
            month=displayTimeCopy.assign(0,spacepos-1);
            displayTimeCopy.erase(0,spacepos);
            spacepos=displayTimeCopy.find_first_of(" ");
            d=displayTimeCopy.assign(0,spacepos-1);
            displayTimeCopy.erase(0,spacepos);
            spacepos=displayTimeCopy.find_first_of(" ");
            h=displayTimeCopy.assign(0,spacepos-1);
            displayTimeCopy.erase(0,spacepos);
            spacepos=displayTimeCopy.find_first_of(" ");
            m=displayTimeCopy.assign(0,spacepos-1);
            displayTimeCopy.erase(0,spacepos);
            spacepos=displayTimeCopy.find_first_of(" ");
            s=displayTimeCopy.assign(0,spacepos-1);
            displayTimeCopy.erase(0,spacepos);
            spacepos=displayTimeCopy.find_first_of(" ");
            month=displayTimeCopy.assign(0,spacepos-1);
            displayTimeCopy.erase(0,spacepos);
            spacepos=displayTimeCopy.find_first_of(" ");
            if(spacepos>10){//no space present; end reached
                    yr=displayTimeCopy;
                    displayTimeCopy="";}
                else{
                    yr=displayTimeCopy.assign(0,spacepos-1);
                    displayTimeCopy.erase(0,spacepos);}
            while(m.length()<2){m.insert(0,"0");}
            while(s.length()<2){m.insert(0,"0");}
                        //order elements: [0]=yr, [1]=mo, [2]=day, [3]=H, [4]=min, [5]=sec, [6]=use am/pm,
                        //[7]=include day of wk in text, [8]=use mo in text
            bool useAmPm=orderToReport[6], useDay = orderToReport[7], useMonth=orderToReport[8];
            int t, mo,t_hr=orderToReport[3],t_min=orderToReport[4],t_sec=orderToReport[5];
            MonthTextToNum(month,mo);
            std::vector<std::string> sequence;
            for(t=0;t<6;++t){
                std::string x;
                sequence.push_back(x);}
            sequence[orderToReport[0]]=yr;
            sequence[orderToReport[1]]=month;
            sequence[orderToReport[2]]=d;
            sequence[orderToReport[3]]=h;
            sequence[orderToReport[4]]=m;
            sequence[orderToReport[5]]=s;
            int hr = stoi(h);//string to int - built-in C++ function
            if(useMonth){ sequence[orderToReport[1]]=month;}
            if(useAmPm){
                if(hr>12) hr -=12;
                std::stringstream sshr;
                sshr<<hr;
                sequence[orderToReport[3]]=sshr.str();
                std::string ampm("pm");
                if(hr<12) ampm="am";
                std::vector<std::string>::iterator it;
                it = sequence.begin()+t_sec+1;
                sequence.insert(it,ampm);
                if(t_hr>t_sec)t_hr++;
                if(t_min>t_sec)t_min++;
                }
            if(useDay){
                std::vector<std::string>::iterator it;
                it = sequence.begin();
                sequence.insert(it,day);
                t_hr++;t_min++;t_sec++;
                }
            std::stringstream dataline;
            for(int t=0;t<sequence.size();++t){
                dataline<<sequence[t];
                if(t==sequence.size()-1) break;
                if((t==t_hr && t_min==t+1) || (t==t_min && t_sec==t+1)){
                        dataline<<colon;}
                    else{dataline<<" ";}
                }//t
            std::string *display=new std::string(dataline.str());
            return display;
            }//localTimeAtMoment()
            

                //order elements: [0]=yr, [1]=mo, [2]=day, [3]=H, [4]=min, [5]=sec, [6]=use am/pm,
                //[7]=include day of wk in text, [8]=use mo in text
                //e.g. output:
                    //012345000 ==> 2013 10/05 15:31:45
                    //012345100 --> 2013 10/05 3:31:45 pm
                    //201345111 --> Sat October 5, 2013 3:31:45 pm
                //order==null: default Sat Oct 5
        
        int daysInMonth(int m, int yr){//4-digit years only; months 1-12 not 0-11
            if(m==2){
                if(yr==2000){return 28;}
                if(MOD(yr,4)==0){return 29;}
                return 28;
                }
            switch (m){
                case 1: return 31;
                case 3: return 31;
                case 5: return 31;
                case 7: return 31;
                case 8: return 31;
                case 10: return 31;
                case 12: return 31;
                default: return 30;
                }
            return 30;
            }//end daysInMonth
        
    std::string* timeStampFilenameFormat(void){//at the moment the function is called
            //20140717_09h00m37s = 7/17/2014 at 9:00:37 am
            //this sorts by runtime
            time_t rightNow = now();
            std::string *timestamp = timeStampFilenameFormat(rightNow);
            return timestamp;
            }
    
    std::string* timeStampFilenameFormat(time_t momentInTime){
        int order[9];order[0]=0;order[1]=1;order[2]=2;order[3]=3;order[4]=4;
        order[5]=5;order[6]=order[7]=order[8]=0;//012345000
        std::string space=" ",slash="/",colon=":",day,ymd,mo,hms,sec,min,daddr_t;
        std::string *timestamp = localTimeAtMoment(momentInTime,order);
        std::string suffixBase = *timestamp;
        size_t spacepos=suffixBase.find_first_of(" ");
        day=suffixBase.assign(0,spacepos-1);
        suffixBase.erase(0,spacepos);
        spacepos=suffixBase.find_first_of(" ");
        ymd=suffixBase.assign(0,spacepos-1);
        suffixBase.erase(0,spacepos);
        spacepos=suffixBase.find_first_of(" ");
        mo=suffixBase.assign(0,spacepos-1);
        suffixBase.erase(0,spacepos);
        spacepos=suffixBase.find_first_of(" ");
        day=suffixBase.assign(0,spacepos-1);
        suffixBase.erase(0,spacepos);
        while(mo.length()<2){mo.insert(0,"0");}
        while(day.length()<2){day.insert(0,"0");}
        suffixBase.append("s");
        size_t colonpos = suffixBase.find_last_of(colon);
        colonpos= suffixBase.find_last_of(colon);
        sec=suffixBase.assign(colonpos,suffixBase.length()-1);
        suffixBase.erase(colonpos,suffixBase.length()-1);
        suffixBase.append("m");
        colonpos= suffixBase.find_last_of(colon);
        min=suffixBase.assign(colonpos,suffixBase.length()-1);
        suffixBase.erase(colonpos,suffixBase.length()-1);
        hms= suffixBase; hms.append("h");
        hms.append(min); hms.append(sec);
        ymd.append(mo);ymd.append(day);
        timestamp->clear();
        timestamp->append(ymd);timestamp->append("_");timestamp->append(hms);
        
        return timestamp;
        }//aTime::timeStampFilenameFormat(time_t )

    
    
        std::string* timeStampFilenamePrefix(void){//for use as a filename prefix
            std::string *filenamebase = timeStampFilenameFormat();
            filenamebase->append("_");
            return filenamebase;}

    std::string* timeStampFilenamePrefix(time_t momentInTime){
            std::string *filenamebase = timeStampFilenameFormat(momentInTime);
            filenamebase->append("_");
            return filenamebase;}

    std::string* timeStampFilenameSuffix(void){//for use as a filename suffix
            std::string *filenamebase = timeStampFilenameFormat();
            std::stringstream s;
            s<<"_"<<*filenamebase;
            *filenamebase=s.str();
            return filenamebase;}

    std::string* timeStampFilenameSuffix(time_t momentInTime){
            std::string *filenamebase = timeStampFilenameFormat(momentInTime);
            std::stringstream s;
            s<<"_"<<*filenamebase;
            *filenamebase=s.str();
            return filenamebase;}

            
    std::string* timeString(time_t momentInTime){
        //string in day of wk/mo/day/H:M:S/yr
        int order[9];order[0]=0;order[1]=1;order[2]=2;order[3]=3;order[4]=4;
        order[5]=5;order[6]=order[7]=order[8]=0;//012345000
        std::string space=" ",slash="/",colon=":",ymd,mo,day,sec,min,hms;
        std::string *timestamp = localTimeAtMoment(momentInTime,order);
        std::string suffixBase = *timestamp;
        size_t spacepos=suffixBase.find_first_of(" ");
        ymd=suffixBase.assign(0,spacepos-1);
        suffixBase.erase(0,spacepos);
        spacepos=suffixBase.find_first_of(" ");
        mo=suffixBase.assign(0,spacepos-1);
        suffixBase.erase(0,spacepos);
        spacepos=suffixBase.find_first_of(" ");
        day=suffixBase.assign(0,spacepos-1);
        suffixBase.erase(0,spacepos);
        while(mo.length()<2){mo.insert(0,"0");}
        while(day.length()<2){day.insert(0,"0");}
        suffixBase.append("s");
        size_t colonpos = suffixBase.find_last_of(colon);
        colonpos= suffixBase.find_last_of(colon);
        sec=suffixBase.assign(colonpos,suffixBase.length()-1);
        suffixBase.erase(colonpos,suffixBase.length()-1);
        suffixBase.append("m");
        colonpos = suffixBase.find_last_of(colon);
        colonpos= suffixBase.find_last_of(colon);
        min=suffixBase.assign(colonpos,suffixBase.length()-1);
        suffixBase.erase(colonpos,suffixBase.length()-1);
        hms= suffixBase; hms.append("h");
        hms.append(min); hms.append(sec);
        ymd.append(mo);ymd.append(day);
        timestamp->clear();
        timestamp->append(ymd);timestamp->append(" ");timestamp->append(hms);
        
        return timestamp;
        }//timeString(time_t )

    long double totalElapsed(void){//firstTime_ to now, in seconds
        time_t now;
        time(&now);
        return difftime(now,firstTime_);}

    friend int operator==(aTime& t1, aTime& t2);
        //true if all data equal
    
    friend int operator!=(aTime& t1, aTime& t2);
        
        };//end class aTime

int operator==(aTime& t1, aTime& t2){
                return ((t1.firstTime_==t2.firstTime_)&&
                        (t1.mostRecentTime_==t2.mostRecentTime_)&&
                        (t1.elapsedSinceLast_==t2.elapsedSinceLast_) );}
        
int operator!=(aTime& t1, aTime& t2){
                return (! (t1==t2) );}


template <class T> void indicesAtMax(T *val[], int items, T* max, int* x, int& numMaxima){
    numMaxima=1;
    max=val[0]; x[0]=0;
    for(int i=1;i<items;++i){
        if(val[i]>max){
                numMaxima=1;
                max=val[i]; x[numMaxima-1]=i;}
            else if(val[i]==max){
                numMaxima++;
                x[numMaxima+1]=i;}
            else{}
        }
    }//indicesAtMax

template <class T> void indicesAtMax(std::vector<T>& val, int items, std::vector<T>& max, std::vector<int>& x, int& numMaxima){
    numMaxima=1;
    max=val[0]; x[0]=0;
    for(int i=1;i<items;++i){
        if(val[i]>max){
                numMaxima=1;
                max=val[i]; x[numMaxima-1]=i;}
            else if(val[i]==max){
                numMaxima++;
                x[numMaxima+1]=i;}
            else{}
        }
    }//indicesAtMax


template <class T> void indicesAtMin(T *val[], int items, T* min, int* x, int& numMinima){
    numMinima=1;
    min=val[0]; x[0]=0;
    for(int i=1;i<items;++i){
        if(val[i]<min){
                numMinima=1;
                min=val[i]; x[numMinima-1]=i;}
            else if(val[i]==min){
                numMinima++;
                x[numMinima-1]=i;}
            else{}
        }
    }//indexAtMin


template <class T> void indicesAtMin(std::vector<T>& val, int items, std::vector<T>& min, std::vector<int>& x, int& numMinima){
    numMinima=1;
    min=val[0]; x[0]=0;
    for(int i=1;i<items;++i){
        if(val[i]<min){
                numMinima=1;
                min=val[i]; x[numMinima-1]=i;}
            else if(val[i]==min){
                numMinima++;
                x[numMinima-1]=i;}
            else{}
        }
    }//indexAtMin



template <class T> void CoordinatesAtMax(T* val, T* xcoords, T* ycoords, int items, T& max, T** xy, int& numMaxima){
    long double p,q;//for debugging
    max=val[0]; p=xy[0][0]=xcoords[0]; q=xy[0][1]=ycoords[0];
    numMaxima=1;
    for(int i=1;i<items;++i){
        if(val[i]>max){
                numMaxima=1;
                max=val[i]; p=xy[0][0]=xcoords[i]; q=xy[0][1]=ycoords[i];}
            else if(val[i]==max){
                numMaxima++;
                p=xy[numMaxima-1][0]=xcoords[i]; q=xy[numMaxima-1][1]=ycoords[i];}
            else{}
        }
    }//CoordinatesAtMax

template <class T> void CoordinatesAtMax(std::vector<T>& val, std::vector<T>& xcoords,
                                         std::vector<T>& ycoords, int items, T& max,
                                         std::vector<std::vector<T> >& xy, int& numMaxima){
long double p,q;//for debugging
    max=val[0]; p=xy[0][0]=xcoords[0]; q=xy[0][1]=ycoords[0];
    numMaxima=1;
    for(int i=1;i<items;++i){
        if(val[i]>max){
                numMaxima=1;
                max=val[i]; p=xy[0][0]=xcoords[i]; q=xy[0][1]=ycoords[i];}
            else if(val[i]==max){
                numMaxima++;
                p=xy[numMaxima-1][0]=xcoords[i]; q=xy[numMaxima-1][1]=ycoords[i];}
            else{}
        }
    }//CoordinatesAtMax



template <class T> void CoordinatesAtMin(T* val, T* xcoords, T* ycoords, int items, T& min, T** xy, int& numMinima){
    min=val[0]; xy[0][0]=xcoords[0]; xy[0][1]=ycoords[0];
    numMinima=1;
    for(int i=1;i<items;++i){
        if(val[i]<min){
                numMinima=1;
                min=val[i]; xy[0][0]=xcoords[i]; xy[0][1]=ycoords[i];}
            else if(val[i]==min){
                numMinima++;
                xy[numMinima-1][0]=xcoords[i]; xy[numMinima-1][1]=ycoords[i];}
            else{}
        }
    }//CoordinatesAtMin

template <class T> void CoordinatesAtMin(std::vector<T>& val, std::vector<T>& xcoords,
                                         std::vector<T>& ycoords, int items, T& min,
                                         std::vector<std::vector<T> >& xy, int& numMinima){
    min=val[0]; xy[0][0]=xcoords[0]; xy[0][1]=ycoords[0];
    numMinima=1;
    for(int i=1;i<items;++i){
        if(val[i]<min){
                numMinima=1;
                min=val[i]; xy[0][0]=xcoords[i]; xy[0][1]=ycoords[i];}
            else if(val[i]==min){
                numMinima++;
                xy[numMinima-1][0]=xcoords[i]; xy[numMinima-1][1]=ycoords[i];}
            else{}
        }
    }//CoordinatesAtMin




template <class T> void IndicesAtMax(T* val, int* xcoords, int* ycoords, int items, T& max, int** xy, int numMaxima){
    numMaxima=1;
    max=val[0]; xy[0][0]=xcoords[0]; xy[0][1]=ycoords[0];
    for(int i=1;i<items;++i){
        if(val[i]>max){
                numMaxima=1;
                max=val[i]; xy[0][0]=xcoords[i]; xy[0][1]=ycoords[i];}
            else if(val[i]==max){
                numMaxima++;
                xy[numMaxima-1][0]=xcoords[i]; xy[numMaxima-1][1]=ycoords[i];}
            else{}
        }
    }//IndicesAtMax

template <class T> void IndicesAtMax(std::vector<T>& val, std::vector<int>& xcoords,
                                     std::vector<int>& ycoords, int items, T& max,
                                     std::vector<std::vector<int> >& xy, int numMaxima){
    numMaxima=1;
    max=val[0]; xy[0][0]=xcoords[0]; xy[0][1]=ycoords[0];
    for(int i=1;i<items;++i){
        if(val[i]>max){
                numMaxima=1;
                max=val[i]; xy[0][0]=xcoords[i]; xy[0][1]=ycoords[i];}
            else if(val[i]==max){
                numMaxima++;
                xy[numMaxima-1][0]=xcoords[i]; xy[numMaxima-1][1]=ycoords[i];}
            else{}
        }
    }//IndicesAtMax



template <class T> void IndicesAtMin(T* val, int* xcoords, int* ycoords, int items, T& min, int** xy, int numMinima){
    numMinima=1;
    min=val[0]; xy[0][0]=xcoords[0]; xy[0][1]=ycoords[0];
    for(int i=1;i<items;++i){
        if(val[i]<min){
                numMinima=1;
                min=val[i]; xy[0][0]=xcoords[i]; xy[0][1]=ycoords[i];}
            else if(val[i]==min){
                numMinima++;
                xy[numMinima-1][0]=xcoords[i]; xy[numMinima-1][1]=ycoords[i];}
            else{}
        }
    }//IndicesAtMin

template <class T> void IndicesAtMin(std::vector<T>& val, std::vector<int>& xcoords,
                                     std::vector<int>& ycoords, int items, T& min,
                                     std::vector<std::vector<int> >& xy, int numMinima){
    numMinima=1;
    min=val[0]; xy[0][0]=xcoords[0]; xy[0][1]=ycoords[0];
    for(int i=1;i<items;++i){
        if(val[i]<min){
                numMinima=1;
                min=val[i]; xy[0][0]=xcoords[i]; xy[0][1]=ycoords[i];}
            else if(val[i]==min){
                numMinima++;
                xy[numMinima-1][0]=xcoords[i]; xy[numMinima-1][1]=ycoords[i];}
            else{}
        }
    }//IndicesAtMin



bool stringVectorIncludesItem(std::vector<std::string>& strvector, std::string& str, long& pos){
    pos=0;
    long items=strvector.size();
    for(long i=0;i<items;++i){
        if(strvector[i]==str){
            pos=i; return true;}}
    return false;
    }//stringVectorIncludesItem


template <class T> void CullDuplicateCoordinates2D(T** coords, int& items){
if(items<2) return;
int i=items-1;
while(i>0){
    for(int j=i-1;j>=0;j--){//from position i, go backwards down the list looking for a duplicate
        //if the ith matches the jth
        if(coords[i][0]==coords[j][0] && coords[i][1] == coords[j][1]){
            //remove the ith coordinate by shifting subsequent items one slot earlier
            if(i<items-1){//if we're removing any but the last item
                for(int k=i;k<items-1;++k){
                    coords[k][0]=coords[k+1][0];
                    coords[k][1]=coords[k+1][1];
                    }
                break;
                }
            items--;
            }
        }//j
    i--;
    }//i
}//CullDuplicateCoordinates2D()


template <class T> void CullDuplicateCoordinates2D(std::vector<std::vector<T> >& coords, int& items){
    if(items<2) return;
    int i=items-1;
    while(i>0){
        for(int j=i-1;j>=0;j--){//from position i, go backwards down the list looking for a duplicate
            //if the ith matches the jth
            if(coords[i][0]==coords[j][0] && coords[i][1] == coords[j][1]){
                //remove the ith coordinate by shifting subsequent items one slot earlier
                if(i<items-1){//if we're removing any but the last item
                    for(int k=i;k<items-1;++k){
                        coords[k][0]=coords[k+1][0];
                        coords[k][1]=coords[k+1][1];
                        }
                    break;
                    }
                items--;
                }
            }//j
        i--;
        }//i
    }//CullDuplicateCoordinates2D()




template <class T, class U> void CullDuplicateCoordinates2D(T** coords, U* pairedList, int& items){
    if(items<2) return;
    int i=items-1;
    while(i>0){
        for(int j=i-1;j>=0;j--){//from position i, go backwards down the list looking for a duplicate
            //if the ith matches the jth
            if(coords[i][0]==coords[j][0] && coords[i][1] == coords[j][1]){
                //remove the ith coordinate by shifting subsequent items one slot earlier
                if(i<items-1){//if we're removing any but the last item
                    for(int k=i;k<items-1;++k){
                        coords[k][0]=coords[k+1][0];
                        coords[k][1]=coords[k+1][1];
                        pairedList[k]=pairedList[k+1];
                        }
                    items--;
                    break;
                    }
                }
            }//j
        i--;
        }//i
    }//CullDuplicateCoordinates2D()


template <class T, class U> void CullDuplicateCoordinates2D(std::vector<std::vector<T> >& coords,
                                                            std::vector<U>& pairedList, int& items){
    if(items<2) return;
    int i=items-1;
    while(i>0){
        for(int j=i-1;j>=0;j--){//from position i, go backwards down the list looking for a duplicate
            //if the ith matches the jth
            if(coords[i][0]==coords[j][0] && coords[i][1] == coords[j][1]){
                //remove the ith coordinate by shifting subsequent items one slot earlier
                if(i<items-1){//if we're removing any but the last item
                    for(int k=i;k<items-1;++k){
                        coords[k][0]=coords[k+1][0];
                        coords[k][1]=coords[k+1][1];
                        pairedList[k]=pairedList[k+1];
                        }
                    items--;
                    break;
                    }
                }
            }//j
        i--;
        }//i
    }//CullDuplicateCoordinates2D()



template <class T> inline int countSetBits(T n) {
	int count=0;
    while( n!= 0){
        n &= (n-1); ++count;}
    return count;}

//Hamming distance counts the number of bits that differ between two binary numbers
//T should be an unsigned integer type
template <class T> inline int HammingDistance(T m, T n) {
	T differentBits = m^n;
    return countSetBits(differentBits);}

enum typeOfModelToRun {dosageOnly=0,tfProductOnly=1,cisOnly=2,tfOnly=3,allSites=4};


class genotypeSettings{
    public:
    uint64_t dosageVal0_,dosageVal1_,tfVal0_,tfVal1_,cisVal0_,cisVal1_;
    public:
    genotypeSettings(void){
        dosageVal0_=dosageVal1_=tfVal0_=tfVal1_=cisVal0_=cisVal1_=0;}
    genotypeSettings(uint64_t dosageVal0, uint64_t dosageVal1, uint64_t tfVal0, uint64_t tfVal1,
                uint64_t cisVal0, uint64_t cisVal1):
                    dosageVal0_(dosageVal0),dosageVal1_(dosageVal1),tfVal0_(tfVal0),tfVal1_(tfVal1),
                    cisVal0_(cisVal0),cisVal1_(cisVal1){}
    genotypeSettings(const genotypeSettings& gs){
        dosageVal0_=gs.dosageVal0_; dosageVal1_=gs.dosageVal1_;
        tfVal0_=gs.tfVal0_; tfVal1_=gs.tfVal1_;
        cisVal0_=gs.cisVal0_; cisVal1_=gs.cisVal1_;}
    ~genotypeSettings(void){}
    genotypeSettings& operator=(const genotypeSettings& gs){
        dosageVal0_=gs.dosageVal0_; dosageVal1_=gs.dosageVal1_;
        tfVal0_=gs.tfVal0_; tfVal1_=gs.tfVal1_;
        cisVal0_=gs.cisVal0_; cisVal1_=gs.cisVal1_;
        return *this;}
        
    };//genotypeSettings




class simulationSettings{
    public:
    int bitstringLen_;
    long double NtfsatPerAllele_,deltaG1dosage_,deltaG1_;
    long double minExpression_,maxExpression_;
    long double Popt_,omega_;
    bool splitSinglePoptRun_;
    uint64_t startingTF0val_, endTF0val_;
    public:
    simulationSettings(void):bitstringLen_(0),splitSinglePoptRun_(false),startingTF0val_(0),endTF0val_(0){
        NtfsatPerAllele_=deltaG1dosage_=deltaG1_=minExpression_=maxExpression_=Popt_=omega_=zero;}
    simulationSettings(int bitstringLen, long double NtfsatPerAllele, long double deltaG1dosage, long double deltaG1,
                        long double minExpression, long double maxExpression, long double Popt, long double omega):
            bitstringLen_(bitstringLen),NtfsatPerAllele_(NtfsatPerAllele),deltaG1dosage_(deltaG1dosage),deltaG1_(deltaG1),
            minExpression_(minExpression),maxExpression_(maxExpression),Popt_(Popt),omega_(omega),
            splitSinglePoptRun_(false),startingTF0val_(0),endTF0val_(0){}
    simulationSettings(int bitstringLen, long double NtfsatPerAllele, long double deltaG1dosage, long double deltaG1,
                        long double minExpression, long double maxExpression, long double Popt, long double omega,
                        bool splitSinglePoptRun,uint64_t startingTF0val, uint64_t endTF0val):
            bitstringLen_(bitstringLen),NtfsatPerAllele_(NtfsatPerAllele),deltaG1dosage_(deltaG1dosage),deltaG1_(deltaG1),
            minExpression_(minExpression),maxExpression_(maxExpression),Popt_(Popt),omega_(omega),
            splitSinglePoptRun_(splitSinglePoptRun),startingTF0val_(startingTF0val),endTF0val_(endTF0val){}
    simulationSettings(const simulationSettings& ss){
        bitstringLen_=ss.bitstringLen_;
        NtfsatPerAllele_=ss.NtfsatPerAllele_; deltaG1dosage_=ss.deltaG1dosage_; deltaG1_=ss.deltaG1_;
        minExpression_=ss.minExpression_; maxExpression_=ss.maxExpression_;
        Popt_=ss.Popt_; omega_=ss.omega_;
        splitSinglePoptRun_=ss.splitSinglePoptRun_; startingTF0val_=ss.startingTF0val_; endTF0val_=ss.endTF0val_;
        }
    ~simulationSettings(void){}
    simulationSettings& operator=(const simulationSettings& ss){
        bitstringLen_=ss.bitstringLen_;
        NtfsatPerAllele_=ss.NtfsatPerAllele_; deltaG1dosage_=ss.deltaG1dosage_; deltaG1_=ss.deltaG1_;
        minExpression_=ss.minExpression_; maxExpression_=ss.maxExpression_;
        Popt_=ss.Popt_; omega_=ss.omega_;
        splitSinglePoptRun_=ss.splitSinglePoptRun_; startingTF0val_=ss.startingTF0val_; endTF0val_=ss.endTF0val_;
        return *this;
        }

    };//simulationSettings





class SimplestRegPathIndividual;
class SimplestRegPathIndividual{
  public:
		uint64_t TFdosage_[2];
		uint64_t TFproduct_[2];
		uint64_t cis_[2];
		int mTFdosage_[2];//mismatches
		int mTF01cis01_[2][2];
		std::string mismatchGtype_;
		long double phenotype_;
		long double fitness_;
		bool phenotypeCalculated_;
		bool fitnessCalculated_;
		bool useMismatchesToCalculatePhenotype_;
  public:
    SimplestRegPathIndividual(void):mismatchGtype_(""){
		TFdosage_[0]=TFdosage_[1]=TFproduct_[0]=TFproduct_[1]=cis_[0]=cis_[1]=0;
		mTFdosage_[0]=mTFdosage_[1]=0;
		mTF01cis01_[0][0]=mTF01cis01_[0][1]=mTF01cis01_[1][0]=mTF01cis01_[1][1]=0;
		phenotype_=fitness_=-one;
		
		useMismatchesToCalculatePhenotype_=true;
		phenotypeCalculated_=fitnessCalculated_=false;
		}
	
	SimplestRegPathIndividual(bool useMismatchesToCalculatePhenotype):mismatchGtype_(""){
		TFdosage_[0]=TFdosage_[1]=TFproduct_[0]=TFproduct_[1]=cis_[0]=cis_[1]=0;
		mTFdosage_[0]=mTFdosage_[1]=0;
		mTF01cis01_[0][0]=mTF01cis01_[0][1]=mTF01cis01_[1][0]=mTF01cis01_[1][1]=0;
		phenotype_=fitness_=-one;
		
		useMismatchesToCalculatePhenotype_=useMismatchesToCalculatePhenotype;
		phenotypeCalculated_=fitnessCalculated_=false;
		}
	
	SimplestRegPathIndividual(const SimplestRegPathIndividual& it){
			*this=it;}
	~SimplestRegPathIndividual(void){}

	SimplestRegPathIndividual& operator=(const SimplestRegPathIndividual& it){
		TFdosage_[0]=it.TFdosage_[0]; TFdosage_[1]=it.TFdosage_[1];
		TFproduct_[0]=it.TFproduct_[0]; TFproduct_[1]=it.TFproduct_[1];
		cis_[0]=it.cis_[0]; cis_[1]=it.cis_[1];
		mTFdosage_[0]=it.mTFdosage_[0];
		mTFdosage_[1]=it.mTFdosage_[1];
		mTF01cis01_[0][0]=it.mTF01cis01_[0][0];
		mTF01cis01_[0][1]=it.mTF01cis01_[0][1];
		mTF01cis01_[1][0]=it.mTF01cis01_[1][0];
		mTF01cis01_[1][1]=it.mTF01cis01_[1][1];
		mismatchGtype_=it.mismatchGtype_;
		phenotype_=it.phenotype_; fitness_=it.fitness_;
		useMismatchesToCalculatePhenotype_=it.useMismatchesToCalculatePhenotype_;
		phenotypeCalculated_=it.phenotypeCalculated_;
		fitnessCalculated_=it.fitnessCalculated_;
		return *this;
		}//operator=


std::string binaryFromUnsignedLongLong(unsigned long long val, bool fullString, int prependZerosToLen){
    int i,last1=0;
    size_t length=CHAR_BIT*sizeof(unsigned long long);
    std::string bitstring(length,'0');
    unsigned int lastbitMask=1;
    if(val==16){
        lastbitMask+=0;}
    for(i=0;i<length;++i){
        bitstring[length-1-i] ='0'+char(val & lastbitMask);
        if(bitstring[length-1-i]=='1') {last1=i+1;}
        val = val>>1;
        }
    if(fullString) return bitstring;
    bitstring.erase(bitstring.begin(),bitstring.end()-last1);
    std::string zeros(prependZerosToLen-last1,'0');
    zeros+=bitstring;
    return zeros;
    }//binaryFromUnsignedLongLong


	std::string gtypeString(int bitstringLen){
		std::string gstr="{{";
		gstr+=binaryFromUnsignedLongLong(TFdosage_[0],false,bitstringLen)+".";
		gstr+=binaryFromUnsignedLongLong(TFproduct_[0],false,bitstringLen)+" || ";
		gstr+=binaryFromUnsignedLongLong(TFdosage_[1],false,bitstringLen)+".";
		gstr+=binaryFromUnsignedLongLong(TFproduct_[1],false,bitstringLen)+"},{";
		gstr+=binaryFromUnsignedLongLong(cis_[0],false,bitstringLen)+" || ";
		gstr+=binaryFromUnsignedLongLong(cis_[1],false,bitstringLen)+"}}";
		return gstr;
		}//gtypeString


	std::string mismatchStringMathematicaFormat(void){
        if(mismatchGtype_[0]=='{'){//it's been calculated
            return mismatchGtype_;}
        std::string gstr="{{";
		uint64_t tacitTF=0;
		if(useMismatchesToCalculatePhenotype_){
				gstr+='0'+char(mTFdosage_[0]);
				gstr+=",";
				gstr+='0'+char(mTFdosage_[1]);
				gstr+="},{{";
				gstr+='0'+char(mTF01cis01_[0][0]);
				gstr+=",";
				gstr+='0'+char(mTF01cis01_[1][0]);
				gstr+="},{";
				gstr+='0'+char(mTF01cis01_[0][1]);
				gstr+=",";
				gstr+='0'+char(mTF01cis01_[1][1]);
				}
			else{
				gstr+='0'+char(HammingDistance(tacitTF,TFdosage_[0]));
				gstr+=",";
				gstr+='0'+char(HammingDistance(tacitTF,TFdosage_[1]));
				gstr+="},{{";
				gstr+='0'+char(HammingDistance(TFproduct_[0],cis_[0]));
				gstr+=",";
				gstr+='0'+char(HammingDistance(TFproduct_[1],cis_[0]));
				gstr+="},{";
				gstr+='0'+char(HammingDistance(TFproduct_[0],cis_[1]));
				gstr+=",";
				gstr+='0'+char(HammingDistance(TFproduct_[1],cis_[1]));
				}
		gstr+="}}}";
		mismatchGtype_=gstr;
		return gstr;
		}//mismatchStringMathematicaFormat
	
	
	long double CalculateFitness(long double Popt,long double omega){
		if(fitnessCalculated_){ return fitness_;}
		fitness_ = exp(-(phenotype_-Popt)*(phenotype_-Popt)/(omega*omega));
		fitnessCalculated_=true;
		return fitness_;
		}

	void CalculateMinMaxExpression(int bitstringLen,long double NtfsatPerAllele, long double deltaG1dosage, long double deltaG1, long double& minExpression, long double& maxExpression, typeOfModelToRun model=allSites){
		long double mDosage0,mDosage1,mTF0cis0,mTF0cis1,mTF1cis0,mTF1cis1;
		mDosage0=mDosage1=mTF0cis0=mTF0cis1=mTF1cis0=mTF1cis1=one;
        if(model==cisOnly || model==tfProductOnly){//dosage at max expression
                mDosage0=mDosage1=zero;}
            else if(model==dosageOnly){//max TF-cis binding
                mTF0cis0=mTF0cis1=mTF1cis0=mTF1cis1=zero;}
            else{}
		long double alphaDose01 = one+NtfsatPerAllele*exp(mDosage0*deltaG1dosage);
		long double alphaDose10 = one+NtfsatPerAllele*exp(mDosage1*deltaG1dosage);
		long double thetaDosage0 = NtfsatPerAllele/(NtfsatPerAllele + alphaDose10*exp(-mDosage0*deltaG1dosage));
		long double thetaDosage1 = NtfsatPerAllele/(NtfsatPerAllele + alphaDose01*exp(-mDosage1*deltaG1dosage));
		long double Ntf0=thetaDosage0*NtfsatPerAllele;
		long double Ntf1=thetaDosage1*NtfsatPerAllele;
		long double alpha00 = one+Ntf0*exp(mTF0cis0*deltaG1);
		long double alpha10 = one+Ntf1*exp(mTF1cis0*deltaG1);
		long double alpha01 = one+Ntf0*exp(mTF0cis1*deltaG1);
		long double alpha11 = one+Ntf1*exp(mTF1cis1*deltaG1);
		long double theta00 = Ntf0/(Ntf0 + alpha10*exp(-mTF0cis0*deltaG1));
		long double theta10 = Ntf1/(Ntf1 + alpha00*exp(-mTF1cis0*deltaG1));
		long double theta01 = Ntf0/(Ntf0 + alpha11*exp(-mTF0cis1*deltaG1));
		long double theta11 = Ntf1/(Ntf1 + alpha01*exp(-mTF1cis1*deltaG1));
		minExpression = (theta00+theta10+theta01+theta11)/two;

		mDosage0=mDosage1=mTF0cis0=mTF0cis1=mTF1cis0=mTF1cis1=zero;
		alphaDose01 = one+NtfsatPerAllele*exp(mDosage0*deltaG1dosage);
		alphaDose10 = one+NtfsatPerAllele*exp(mDosage1*deltaG1dosage);
		thetaDosage0 = NtfsatPerAllele/(NtfsatPerAllele + alphaDose10*exp(-mDosage0*deltaG1dosage));
		thetaDosage1 = NtfsatPerAllele/(NtfsatPerAllele + alphaDose01*exp(-mDosage1*deltaG1dosage));
		Ntf0=thetaDosage0*NtfsatPerAllele;
		Ntf1=thetaDosage1*NtfsatPerAllele;
		alpha00 = one+Ntf0*exp(mTF0cis0*deltaG1);
		alpha10 = one+Ntf1*exp(mTF1cis0*deltaG1);
		alpha01 = one+Ntf0*exp(mTF0cis1*deltaG1);
		alpha11 = one+Ntf1*exp(mTF1cis1*deltaG1);
		theta00 = Ntf0/(Ntf0 + alpha10*exp(-mTF0cis0*deltaG1));
		theta10 = Ntf1/(Ntf1 + alpha00*exp(-mTF1cis0*deltaG1));
		theta01 = Ntf0/(Ntf0 + alpha11*exp(-mTF0cis1*deltaG1));
		theta11 = Ntf1/(Ntf1 + alpha01*exp(-mTF1cis1*deltaG1));
		maxExpression = (theta00+theta10+theta01+theta11)/two;
		}//CalculateMinMaxExpression


	long double CalculatePhenotype(const simulationSettings& ss){
        return CalculatePhenotype(ss.bitstringLen_, ss.NtfsatPerAllele_, ss.deltaG1dosage_,
                    ss.deltaG1_, ss.minExpression_, ss.maxExpression_);
        }//CalculatePhenotype()

	long double CalculatePhenotype(int bitstringLen, long double NtfsatPerAllele, long double deltaG1dosage, long double deltaG1, long double minExpression, long double maxExpression){
		if(phenotypeCalculated_){return phenotype_;}
		long double mDosage0,mDosage1,mTF0cis0,mTF0cis1,mTF1cis0,mTF1cis1;
		uint64_t tacitTF=(uint64_t) 0;
		if(useMismatchesToCalculatePhenotype_){
				mDosage0=(long double)mTFdosage_[0]/(long double)bitstringLen;
				mDosage1=(long double)mTFdosage_[1]/(long double)bitstringLen;
				mTF0cis0=(long double)mTF01cis01_[0][0]/(long double)bitstringLen;
				mTF0cis1=(long double)mTF01cis01_[0][1]/(long double)bitstringLen;
				mTF1cis0=(long double)mTF01cis01_[1][0]/(long double)bitstringLen;
				mTF1cis1=(long double)mTF01cis01_[1][1]/(long double)bitstringLen;
				}
			else{
				mDosage0=(long double)HammingDistance(tacitTF,TFdosage_[0])/(long double)bitstringLen;
				mDosage1=(long double)HammingDistance(tacitTF,TFdosage_[1])/(long double)bitstringLen;
				mTF0cis0=(long double)HammingDistance(TFproduct_[0],cis_[0])/(long double)bitstringLen;
				mTF0cis1=(long double)HammingDistance(TFproduct_[0],cis_[1])/(long double)bitstringLen;
				mTF1cis0=(long double)HammingDistance(TFproduct_[1],cis_[0])/(long double)bitstringLen;
				mTF1cis1=(long double)HammingDistance(TFproduct_[1],cis_[1])/(long double)bitstringLen;
				}
		long double alphaDose01 = one+NtfsatPerAllele*exp(mDosage0*deltaG1dosage);
		long double alphaDose10 = one+NtfsatPerAllele*exp(mDosage1*deltaG1dosage);
		long double thetaDosage0 = NtfsatPerAllele/(NtfsatPerAllele + alphaDose10*exp(-mDosage0*deltaG1dosage));
		long double thetaDosage1 = NtfsatPerAllele/(NtfsatPerAllele + alphaDose01*exp(-mDosage1*deltaG1dosage));
		long double Ntf0=thetaDosage0*NtfsatPerAllele;
		long double Ntf1=thetaDosage1*NtfsatPerAllele;
		long double alpha00 = one+Ntf0*exp(mTF0cis0*deltaG1);
		long double alpha10 = one+Ntf1*exp(mTF1cis0*deltaG1);
		long double alpha01 = one+Ntf0*exp(mTF0cis1*deltaG1);
		long double alpha11 = one+Ntf1*exp(mTF1cis1*deltaG1);
		long double theta00 = Ntf0/(Ntf0 + alpha10*exp(-mTF0cis0*deltaG1));
		long double theta10 = Ntf1/(Ntf1 + alpha00*exp(-mTF1cis0*deltaG1));
		long double theta01 = Ntf0/(Ntf0 + alpha11*exp(-mTF0cis1*deltaG1));
		long double theta11 = Ntf1/(Ntf1 + alpha01*exp(-mTF1cis1*deltaG1));
		long double thetaUnscaled = (theta00+theta10+theta01+theta11)/two;
		long double scaledExpression = (thetaUnscaled-minExpression)/(maxExpression-minExpression);
		phenotype_ = MAX(scaledExpression,zero);
		phenotypeCalculated_=true;
		return phenotype_;
		}//CalculatePhenotype

	long double fitness(void){return fitness_;}
	
	std::string hetType(bool useMismatches=false){
		std::string ht="___";
		if(IsTFdosageHeterozygote(useMismatches) && IsTFprodHeterozygote(useMismatches) && IsCisHeterozygote(useMismatches)){
				ht="dpc";}
			else if(IsTFdosageHeterozygote(useMismatches) && IsTFprodHeterozygote(useMismatches) && !(IsCisHeterozygote(useMismatches))){
				ht="dp_";}
			else if(IsTFdosageHeterozygote(useMismatches) && !(IsTFprodHeterozygote(useMismatches)) && IsCisHeterozygote(useMismatches)){
				ht="d_c";}
			else if(IsTFdosageHeterozygote(useMismatches) && !(IsTFprodHeterozygote(useMismatches)) && !(IsCisHeterozygote(useMismatches))){
				ht="d__";}
			else if(!(IsTFdosageHeterozygote(useMismatches)) && IsTFprodHeterozygote(useMismatches) && IsCisHeterozygote(useMismatches)){
				ht="_pc";}
			else if(!(IsTFdosageHeterozygote(useMismatches)) && IsTFprodHeterozygote(useMismatches) && !(IsCisHeterozygote(useMismatches))){
				ht="_p_";}
			else if(!(IsTFdosageHeterozygote(useMismatches)) && !(IsTFprodHeterozygote(useMismatches)) && IsCisHeterozygote(useMismatches)){
				ht="__c";}
			else{}
		return ht;
		}


	bool IsCisHeterozygote(bool useMismatches=false){
		if(useMismatches || useMismatchesToCalculatePhenotype_){
				return (mTF01cis01_[0][0] != mTF01cis01_[0][1]) || (mTF01cis01_[1][0] != mTF01cis01_[1][1]);}
			else{
				return cis_[0]!=cis_[1];}
			}//IsCisHeterozygote

	bool IsTFheterozygote(bool useMismatches=false){
		if(useMismatches || useMismatchesToCalculatePhenotype_){
				return (IsTFdosageHeterozygote() || IsTFprodHeterozygote());}
			else{
				return (TFdosage_[0]!=TFdosage_[1] || TFproduct_[0]!=TFproduct_[1]);}
		}//IsTFheterozygote

	bool IsTFdosageHeterozygote(bool useMismatches=false){
		if(useMismatches || useMismatchesToCalculatePhenotype_){
				return (mTFdosage_[0]!=mTFdosage_[1]);}
			else{
				return (TFdosage_[0]!=TFdosage_[1]);}
		}//IsTFdosageHeterozygote
		
	bool IsTFprodHeterozygote(bool useMismatches=false){
		if(useMismatches || useMismatchesToCalculatePhenotype_){//mathematica: (mTF1cis1 != mTF2cis1) || ( mTF1cis2 != mTF2cis2)
				return (mTF01cis01_[0][0] != mTF01cis01_[1][0]) || (mTF01cis01_[0][1] != mTF01cis01_[1][1]);}
			else{
			return (TFproduct_[0]!=TFproduct_[1]);}
		}//IsTFprodHeterozygote()

    std::string& mismatchGenotype(void){return mismatchGtype_;}
    std::string mismatchHetType(void){return hetType(true);}

    long double phenotype(void){return phenotype_;}
    void Reset(void){phenotypeCalculated_=fitnessCalculated_=false; phenotype_=fitness_=-one;}

	void SetGenotype(int dosageTF, int prodTF, int cisSite, int mismatchVal){
		if(dosageTF==0 || dosageTF==1){//dosage -- set dosageTF out of 0-1 range to set TF-cis
				mTFdosage_[dosageTF]=mismatchVal;}
			else{
				mTF01cis01_[prodTF][cisSite]=mismatchVal;}
		mismatchStringMathematicaFormat();
		phenotypeCalculated_=fitnessCalculated_=false;
		}//SetGenotype
	
	void SetGenotype(int site, int alleleCopy, uint64_t binaryVal){
		if(site==0){//dosage
				TFdosage_[alleleCopy]=binaryVal;
				}
			else if (site==1){//TF product
				TFproduct_[alleleCopy]=binaryVal;
				}
			else{//cis
				cis_[alleleCopy]=binaryVal;}
		mismatchStringMathematicaFormat();
		phenotypeCalculated_=fitnessCalculated_=false;
		}//SetGenotype
	
    void SetMismatchesUsingBitstrings(void){
        uint64_t tacitTF=(uint64_t) 0;
        mTFdosage_[0]=HammingDistance(tacitTF,TFdosage_[0]);
        mTFdosage_[1]=HammingDistance(tacitTF,TFdosage_[1]);
        mTF01cis01_[0][0]=HammingDistance(TFproduct_[0],cis_[0]);
        mTF01cis01_[0][1]=HammingDistance(TFproduct_[0],cis_[1]);
        mTF01cis01_[1][0]=HammingDistance(TFproduct_[1],cis_[0]);
        mTF01cis01_[1][1]=HammingDistance(TFproduct_[1],cis_[1]);
        }//SetMismatchesUsingBitstrings
        
        
	friend int operator==(SimplestRegPathIndividual& i1, SimplestRegPathIndividual& i2);
	friend int operator!=(SimplestRegPathIndividual& i1, SimplestRegPathIndividual& i2);

	};// SimplestRegPathIndividual


int operator==(SimplestRegPathIndividual& i1, SimplestRegPathIndividual& i2){
		if(i1.useMismatchesToCalculatePhenotype_){
				if(i1.mTFdosage_[0]==i2.mTFdosage_[0] && i1.mTFdosage_[1]==i2.mTFdosage_[1]
					&& i1.mTF01cis01_[0][0]==i2.mTF01cis01_[0][0] && i1.mTF01cis01_[0][1]==i2.mTF01cis01_[0][1]
					&& i1.mTF01cis01_[1][0]==i2.mTF01cis01_[1][0] && i1.mTF01cis01_[1][1]==i2.mTF01cis01_[1][1]){
						return 1;}
				}
			else{
				if(i1.TFdosage_[0]==i2.TFdosage_[0] && i1.TFdosage_[1]==i2.TFdosage_[1]
					&& i1.TFproduct_[0]==i2.TFproduct_[0] && i1.TFproduct_[1]==i2.TFproduct_[1]
					&& i1.cis_[0]==i2.cis_[0] && i1.cis_[1]==i2.cis_[1]){return 1;}
				}
		return 0;}

int operator!=(SimplestRegPathIndividual& i1, SimplestRegPathIndividual& i2){
	return !(i1==i2);}






  int codeToInt(std::string& hetcode){
       int code=0;//"___" default
       if(hetcode=="___"){code=0;}
           else if(hetcode=="__c"){code=1;}
           else if(hetcode=="_p_"){code=2;}
           else if(hetcode=="d__"){code=3;}
           else if(hetcode=="_pc"){code=4;}
           else if(hetcode=="d_c"){code=5;}
           else if(hetcode=="dp_"){code=6;}
           else if(hetcode=="dpc"){code=7;}

           else if(hetcode=="__n"){code=8;}
           else if(hetcode=="_n_"){code=9;}
           else if(hetcode=="n__"){code=10;}
           else if(hetcode=="_nn"){code=11;}
           else if(hetcode=="n_n"){code=12;}
           else if(hetcode=="nn_"){code=13;}
           else if(hetcode=="nnn"){code=14;}

           else if(hetcode=="_nc"){code=15;}
           else if(hetcode=="n_c"){code=16;}
           else if(hetcode=="nnc"){code=17;}

           else if(hetcode=="_pn"){code=18;}
           else if(hetcode=="np_"){code=19;}
           else if(hetcode=="npn"){code=20;}

           else if(hetcode=="d_n"){code=21;}
           else if(hetcode=="dn_"){code=22;}
           else if(hetcode=="dnn"){code=23;}

           else if(hetcode=="npc"){code=24;}
           else if(hetcode=="dnc"){code=25;}
           else if(hetcode=="dpn"){code=26;}

           else{}
       return code;
      }//codeToInt()
   
  std::string intToHetcode(int& hetcodeVal){
      std::string hetcode;// default
      switch(hetcodeVal){
          case 0:
              hetcode="___"; break;
          case 1:
              hetcode="__c"; break;
          case 2:
              hetcode="_p_"; break;
          case 3:
              hetcode="d__"; break;
          case 4:
              hetcode="_pc"; break;
          case 5:
              hetcode="d_c"; break;
          case 6:
              hetcode="dp_"; break;
          case 7:
              hetcode="dpc"; break;
          case 8:
              hetcode="__n"; break;
          case 9:
              hetcode="_n_"; break;
          case 10:
              hetcode="n__"; break;
          case 11:
              hetcode="_nn"; break;
          case 12:
              hetcode="n_n"; break;
          case 13:
              hetcode="nn_"; break;
          case 14:
              hetcode="nnn"; break;
          case 15:
              hetcode="_nc"; break;
          case 16:
              hetcode="n_c"; break;
          case 17:
              hetcode="nnc"; break;
          case 18:
              hetcode="_pn"; break;
          case 19:
              hetcode="np_"; break;
          case 20:
              hetcode="npn"; break;
          case 21:
              hetcode="d_n"; break;
          case 22:
              hetcode="dn_"; break;
          case 23:
              hetcode="dnn"; break;
          case 24:
              hetcode="npc"; break;
          case 25:
              hetcode="dnc"; break;
          case 26:
              hetcode="dpn"; break;
          default: hetcode="error in intToHetcode()"; break;
          }//switch
      return hetcode;
     }//intToHetcode()
  
  
  






class BitstringGenotypeData{
    //this class stores a bitstring, its mismatch combination and its heterozygote code
    //as well as its phenotype and fitness
  public:
    std::string bitstringGtype_;
    std::string mismatchPattern_;
    std::string trueHetcode_;
    std::string mismatchHetcode_;
    long double phenotype_;
    long double fitness_;
    bool phenotypeCalculated_;
    bool fitnessCalculated_;
    
  public:
    BitstringGenotypeData(void){}
    BitstringGenotypeData(const BitstringGenotypeData& bgd){
        operator=(bgd);}
    BitstringGenotypeData(std::string bitstringGtype, std::string mismatchPattern, std::string hetcode, std::string& mismatchHetcode):
             bitstringGtype_(bitstringGtype), mismatchPattern_(mismatchPattern),trueHetcode_(hetcode),mismatchHetcode_(mismatchHetcode),
             phenotype_(zero),fitness_(zero),phenotypeCalculated_(false),fitnessCalculated_(false){}
    BitstringGenotypeData(std::string bitstringGtype, std::string mismatchPattern, std::string hetcode, std::string& mismatchHetcode,
                                        long double phenotype):
             bitstringGtype_(bitstringGtype), mismatchPattern_(mismatchPattern),trueHetcode_(hetcode),mismatchHetcode_(mismatchHetcode),
                phenotype_(phenotype),fitness_(-one),phenotypeCalculated_(true),fitnessCalculated_(false){}
    BitstringGenotypeData(std::string bitstringGtype, std::string mismatchPattern, std::string hetcode, std::string& mismatchHetcode,
                                        long double phenotype, long double fitness):
             bitstringGtype_(bitstringGtype), mismatchPattern_(mismatchPattern),trueHetcode_(hetcode),mismatchHetcode_(mismatchHetcode),
                phenotype_(phenotype),fitness_(fitness),phenotypeCalculated_(true),fitnessCalculated_(true){}
    ~BitstringGenotypeData(void){}
    
    BitstringGenotypeData& operator=(const BitstringGenotypeData& bgd){
        bitstringGtype_=bgd.bitstringGtype_;
        mismatchPattern_=bgd.mismatchPattern_;
        trueHetcode_=bgd.trueHetcode_;
        mismatchHetcode_=bgd.mismatchHetcode_;
        phenotype_=bgd.phenotype_;
        fitness_=bgd.fitness_;
        phenotypeCalculated_=bgd.phenotypeCalculated_;
        fitnessCalculated_=bgd.fitnessCalculated_;
        return *this;}
    
    std::string& genotype(void){return bitstringGtype_;}
    std::string& hetcode(void){return trueHetcode_;}
    std::string& mismatchPattern(void){return mismatchPattern_;}
    std::string& mismatchHetcode(void){return mismatchHetcode_;}
    
    bool sameMismatchPattern(std::string& mismatchPattern){//doesn't check for allele order
        return mismatchPattern_==mismatchPattern;}

    void CollectData(SimplestRegPathIndividual& indiv, int bitstringlen){
        indiv.SetMismatchesUsingBitstrings();
        bitstringGtype_=indiv.gtypeString(bitstringlen);
        mismatchPattern_=indiv.mismatchStringMathematicaFormat();
        trueHetcode_=indiv.hetType(false);
        mismatchHetcode_=indiv.hetType(true);
        if(indiv.phenotypeCalculated_==true){
                phenotype_=indiv.phenotype();phenotypeCalculated_=true;}
            else{phenotype_=zero; phenotypeCalculated_=false;}
        if(indiv.fitnessCalculated_==true){
                fitness_=indiv.fitness();fitnessCalculated_=true;}
            else{fitness_=zero; fitnessCalculated_=false;}
        }//CollectData()
    
    void SetBitstring(std::string& bitstringGtype, std::string& hetcode){
        bitstringGtype_=bitstringGtype; trueHetcode_=hetcode;}
    
    void SetFitness(long double& fitness){
        fitness_=fitness; fitnessCalculated_=true;}
    
    void SetHetcode(std::string& hetcode){
        trueHetcode_=hetcode;}
    
    void SetMismatchPattern(std::string& mismatchPattern){
        mismatchPattern_=mismatchPattern;}
    
    void SetMismatchHetcode(std::string& mismatchHetcode){
        mismatchHetcode_=mismatchHetcode;}
    
    void SetPhenotype(long double& phenotype){
        if(phenotype_ == phenotype){
                phenotypeCalculated_=true;}
            else{
                phenotype_=phenotype; phenotypeCalculated_=true;
                fitness_=zero; fitnessCalculated_=false;
                }
        }//SetPhenotype

    int hetcodeToInt(void){return codeToInt(trueHetcode_);}
    int mismatchHetcodeToInt(void){return codeToInt(mismatchHetcode_);}

    

    void Clear(void){
        bitstringGtype_.clear();mismatchPattern_.clear();
        trueHetcode_.clear();mismatchHetcode_.clear();
        phenotype_=fitness_=zero; phenotypeCalculated_=fitnessCalculated_=false;}
    
    friend int operator==(BitstringGenotypeData& bgd1, BitstringGenotypeData& bgd2);
    friend int operator!=(BitstringGenotypeData& bgd1, BitstringGenotypeData& bgd2);
    };//BitstringGenotypeData

int operator==(BitstringGenotypeData& bgd1, BitstringGenotypeData& bgd2){
    return (bgd1.bitstringGtype_ == bgd2.bitstringGtype_);}

int operator!=(BitstringGenotypeData& bgd1, BitstringGenotypeData& bgd2){
    return !(bgd1 == bgd2);}



class FitnessLandscapeParameters{
  public:
    long double wBarMax_,pBarAtMax_;
    long double wAABB_,wAABb_,wAAbb_,wAaBB_,wAaBb_,wAabb_,waaBB_,waaBb_,waabb_;//genotype-specific fitnesses
    long double pAABB_,pAABb_,pAAbb_,pAaBB_,pAaBb_,pAabb_,paaBB_,paaBb_,paabb_;//genotype-specific phenotypes
    bool wBarMaxCalculated_,pBarAtMaxCalculated_;
    static const long double tol_;
  public:
    FitnessLandscapeParameters(void){
        wBarMax_=pBarAtMax_=-one;
        wAABB_=wAABb_=wAAbb_=wAaBB_=wAaBb_=wAabb_=waaBB_=waaBb_=waabb_=-one;
        pAABB_=pAABb_=pAAbb_=pAaBB_=pAaBb_=pAabb_=paaBB_=paaBb_=paabb_=-one;
        wBarMaxCalculated_=pBarAtMaxCalculated_=false;}
    FitnessLandscapeParameters(long double wBarMax, long double pBarAtMax,
                     long double pAABB,long double pAABb,long double pAAbb,
                     long double pAaBB,long double pAaBb,long double pAabb,
                     long double paaBB,long double paaBb,long double paabb,
                     long double wAABB,long double wAABb,long double wAAbb,
                     long double wAaBB,long double wAaBb,long double wAabb,
                    long double waaBB,long double waaBb,long double waabb):wBarMax_(wBarMax),pBarAtMax_(pBarAtMax),
    pAABB_(pAABB),pAABb_(pAABb),pAAbb_(pAAbb),pAaBB_(pAaBB),pAaBb_(pAaBb),pAabb_(pAabb),paaBB_(paaBB),paaBb_(paaBb),paabb_(paabb),
    wAABB_(wAABB),wAABb_(wAABb),wAAbb_(wAAbb),wAaBB_(wAaBB),wAaBb_(wAaBb),wAabb_(wAabb),waaBB_(waaBB),waaBb_(waaBb),waabb_(waabb){}
        
    FitnessLandscapeParameters(const FitnessLandscapeParameters& flp){
        operator=(flp);}
    ~FitnessLandscapeParameters(void){}
    FitnessLandscapeParameters& operator=(const FitnessLandscapeParameters& flp){
        wBarMax_=flp.wBarMax_;
        pBarAtMax_=flp.pBarAtMax_;
        wAABB_=flp.wAABB_;wAABb_=flp.wAABb_;wAAbb_=flp.wAAbb_;
        wAaBB_=flp.wAaBB_;wAaBb_=flp.wAaBb_;wAabb_=flp.wAabb_;
        waaBB_=flp.waaBB_;waaBb_=flp.waaBb_;waabb_=flp.waabb_;
        pAABB_=flp.pAABB_;pAABb_=flp.pAABb_;pAAbb_=flp.pAAbb_;
        pAaBB_=flp.pAaBB_;pAaBb_=flp.pAaBb_;pAabb_=flp.pAabb_;
        paaBB_=flp.paaBB_;paaBb_=flp.paaBb_;paabb_=flp.paabb_;
        wBarMaxCalculated_=flp.wBarMaxCalculated_;
        pBarAtMaxCalculated_=flp.pBarAtMaxCalculated_;
        return *this;}
    
    void Reset(void){
        wBarMax_=pBarAtMax_=-one;
        wAABB_=wAABb_=wAAbb_=wAaBB_=wAaBb_=wAabb_=waaBB_=waaBb_=waabb_=-one;
        pAABB_=pAABb_=pAAbb_=pAaBB_=pAaBb_=pAabb_=paaBB_=paaBb_=paabb_=-one;
        wBarMaxCalculated_=pBarAtMaxCalculated_=false;}

    
    void SetFitnessParameters(long double wAABB,long double wAABb,long double wAAbb,
                              long double wAaBB,long double wAaBb,long double wAabb,
                              long double waaBB,long double waaBb,long double waabb){
        wAABB_=wAABB;wAABb_=wAABb;wAAbb_=wAAbb;
        wAaBB_=wAaBB;wAaBb_=wAaBb;wAabb_=wAabb;
        waaBB_=waaBB;waaBb_=waaBb;waabb_=waabb;}
    
    void SetpBarAtMax(long double pBarAtMax){
        pBarAtMax_=pBarAtMax;pBarAtMaxCalculated_=true;}
    
    void SetPhenotypicParameters(long double pAABB,long double pAABb,long double pAAbb,
                              long double pAaBB,long double pAaBb,long double pAabb,
                              long double paaBB,long double paaBb,long double paabb){
        pAABB_=pAABB;pAABb_=pAABb;pAAbb_=pAAbb;
        pAaBB_=pAaBB;pAaBb_=pAaBb;pAabb_=pAabb;
        paaBB_=paaBB;paaBb_=paaBb;paabb_=paabb;}
    
    void SetWbarMax(long double wBarMax){
        wBarMax_=wBarMax;wBarMaxCalculated_=true;}

    friend int operator==(const FitnessLandscapeParameters& flp1, const FitnessLandscapeParameters& flp2);
    friend int operator!=(const FitnessLandscapeParameters& flp1, const FitnessLandscapeParameters& flp2);
    };// class FitnessLandscapeParameters

const long double FitnessLandscapeParameters::tol_ = 0.000001;


int operator==(const FitnessLandscapeParameters& flp1, const FitnessLandscapeParameters& flp2){
    return(//flp1.wBarMax_==flp2.wBarMax_ && flp2.pBarAtMax_==flp2.pBarAtMax_ &&
           flp1.wAABB_==flp2.wAABB_ && flp1.wAABb_==flp2.wAABb_ && flp1.wAAbb_==flp2.wAAbb_
           && flp1.wAaBB_==flp2.wAaBB_ && flp1.wAaBb_==flp2.wAaBb_ && flp1.wAabb_==flp2.wAabb_
           && flp1.waaBB_==flp2.waaBB_ && flp1.waaBb_==flp2.waaBb_ && flp1.waabb_==flp2.waabb_

           && flp1.pAABB_==flp2.pAABB_ && flp1.pAABb_==flp2.pAABb_ && flp1.pAAbb_==flp2.pAAbb_
           && flp1.pAaBB_==flp2.pAaBB_ && flp1.pAaBb_==flp2.pAaBb_ && flp1.pAabb_==flp2.pAabb_
           && flp1.paaBB_==flp2.paaBB_ && flp1.paaBb_==flp2.paaBb_ && flp1.paabb_==flp2.paabb_
           && flp1.wBarMaxCalculated_==flp2.wBarMaxCalculated_
           && flp1.pBarAtMaxCalculated_==flp2.pBarAtMaxCalculated_
           );}

int operator!=(const FitnessLandscapeParameters& flp1, const FitnessLandscapeParameters& flp2){
    return !(flp1==flp2);}



class FitnessMaximumBitstringSolution{
    //this class stores identical bitstring-genotype solutions that maximize wBar,
    //and the list of reference genotypes that give this solution
    //If there is more than one solution for a given optimization,
    //then use a list of objects of this class
public:
    BitstringGenotypeData solution_;//stores a heterozygote genotype if the solution is polymorphic
    long double wBarMax_;//error code is -one
    long double phBar_;//equilibrium mean phenotype; error codes are -one
    std::vector<uint64_t> TFdosageAlleles_;
    std::vector<uint64_t> TFprodAlleles_;
    std::vector<uint64_t> cisAlleles_;
    long double phat_;//equilibrium TF allele frequency
    long double qhat_;//equilibrium cis allele frequency
    bool pNeutral_;//TF neutral
    bool qNeutral_;//cis neutral
    std::string mismatchPattern_;//mathematica string format
    std::string trueHetcode_;
    std::string mismatchHetcode_;
    //use these for checking fitness surfaces in Mathematica
    std::vector<FitnessLandscapeParameters> fitnessLandscapesAndGPmaps_;
    std::vector<BitstringGenotypeData> referenceGenotypes_;
    bool headerPrinted_;
    static const long double tol_;

public:
    FitnessMaximumBitstringSolution(void){
        wBarMax_=phBar_=phat_=qhat_= -one;//error code
        pNeutral_=qNeutral_=headerPrinted_=false;
        }
    FitnessMaximumBitstringSolution(const FitnessMaximumBitstringSolution& fd){
        solution_=fd.solution_;
        wBarMax_=fd.wBarMax_;
        phBar_=fd.phBar_;
        TFdosageAlleles_=fd.TFdosageAlleles_;
        TFprodAlleles_=fd.TFprodAlleles_;
        cisAlleles_=fd.cisAlleles_;
        phat_=fd.phat_;
        qhat_=fd.qhat_;
        pNeutral_=fd.pNeutral_;
        qNeutral_=fd.qNeutral_;
        fitnessLandscapesAndGPmaps_=fd.fitnessLandscapesAndGPmaps_;
        referenceGenotypes_=fd.referenceGenotypes_;
        headerPrinted_=fd.headerPrinted_;
        }
    ~FitnessMaximumBitstringSolution(void){
        fitnessLandscapesAndGPmaps_.clear();
        referenceGenotypes_.clear();
        }
    
    FitnessMaximumBitstringSolution& operator=(const FitnessMaximumBitstringSolution& fd){
        solution_=fd.solution_;
        wBarMax_=fd.wBarMax_;
        phBar_=fd.phBar_;
        TFdosageAlleles_=fd.TFdosageAlleles_;
        TFprodAlleles_=fd.TFprodAlleles_;
        cisAlleles_=fd.cisAlleles_;
        phat_=fd.phat_;
        qhat_=fd.qhat_;
        pNeutral_=fd.pNeutral_;
        qNeutral_=fd.qNeutral_;
        fitnessLandscapesAndGPmaps_=fd.fitnessLandscapesAndGPmaps_;
        referenceGenotypes_=fd.referenceGenotypes_;
        headerPrinted_=fd.headerPrinted_;
        return *this;
        }//operator=
    

        
    void AppendReferenceGenotype(BitstringGenotypeData& refBGD){
        //doesn't check whether refBGD is equal or equivalent to an already listed genotype
        //i.e., if AaBb is already listed, then whether it refBGD is AaBb, AabB, aABb or aAbB
        referenceGenotypes_.push_back(refBGD);}
    
    void AppendReferenceGenotype(FitnessMaximumBitstringSolution& alternateSolution){
        //use after testing that the alternateSolution is the same as this solution
        AppendReferenceGenotype(alternateSolution.solution_);}
    
    std::string& hetcode(void){return solution_.hetcode();}
    std::string& mismatchHetcode(void){return solution_.mismatchHetcode();}

    bool isSameSolutionMismatchPattern(std::string& testMismatchPattern){
        return solution_.sameMismatchPattern(testMismatchPattern);}

    bool pNeutral(void) {return pNeutral_;}
    bool qNeutral(void) {return qNeutral_;}

    std::string& solutionGenotype(void){
        return solution_.genotype();}
    
    std::string& solutionHetcode(void){
        return solution_.hetcode();}
    
    int solutionHetcodeToInt(void){
        return solution_.hetcodeToInt();}
    
    std::string& solutionMismatchPattern(void){
        return solution_.mismatchPattern();}
    
    std::string& solutionMismatchHetcode(void){
        return solution_.mismatchHetcode();}
    
    int solutionMismatchHetcodeToInt(void){
        return solution_.mismatchHetcodeToInt();}

    
    void PrintLandscapeAndGPmap(std::ostream& outfile, int mapToPrint){
        PrintLandscapeAndGPmapHeader(outfile);
        std::string tab("\t");
        FitnessLandscapeParameters& s = fitnessLandscapesAndGPmaps_[mapToPrint];
        if(s.wAABB_!=-one){outfile<<ROUND(s.wAABB_,decimalDigitsToRound);};outfile<<tab;
        if(s.wAABb_!=-one){outfile<<ROUND(s.wAABb_,decimalDigitsToRound);};outfile<<tab;
        if(s.wAAbb_!=-one){outfile<<ROUND(s.wAAbb_,decimalDigitsToRound);};outfile<<tab;
        if(s.wAaBB_!=-one){outfile<<ROUND(s.wAaBB_,decimalDigitsToRound);};outfile<<tab;
        if(s.wAaBb_!=-one){outfile<<ROUND(s.wAaBb_,decimalDigitsToRound);};outfile<<tab;
        if(s.wAabb_!=-one){outfile<<ROUND(s.wAabb_,decimalDigitsToRound);};outfile<<tab;
        if(s.waaBB_!=-one){outfile<<ROUND(s.waaBB_,decimalDigitsToRound);};outfile<<tab;
        if(s.waaBb_!=-one){outfile<<ROUND(s.waaBb_,decimalDigitsToRound);};outfile<<tab;
        if(s.waabb_!=-one){outfile<<ROUND(s.waabb_,decimalDigitsToRound);};outfile<<tab;

        if(s.pAABB_!=-one){outfile<<ROUND(s.pAABB_,decimalDigitsToRound);};outfile<<tab;
        if(s.pAABb_!=-one){outfile<<ROUND(s.pAABb_,decimalDigitsToRound);};outfile<<tab;
        if(s.pAAbb_!=-one){outfile<<ROUND(s.pAAbb_,decimalDigitsToRound);};outfile<<tab;
        if(s.pAaBB_!=-one){outfile<<ROUND(s.pAaBB_,decimalDigitsToRound);};outfile<<tab;
        if(s.pAaBb_!=-one){outfile<<ROUND(s.pAaBb_,decimalDigitsToRound);};outfile<<tab;
        if(s.pAabb_!=-one){outfile<<ROUND(s.pAabb_,decimalDigitsToRound);};outfile<<tab;
        if(s.paaBB_!=-one){outfile<<ROUND(s.paaBB_,decimalDigitsToRound);};outfile<<tab;
        if(s.paaBb_!=-one){outfile<<ROUND(s.paaBb_,decimalDigitsToRound);};outfile<<tab;
        if(s.paabb_!=-one){outfile<<ROUND(s.paabb_,decimalDigitsToRound);};outfile<<tab;
        }//PrintLandscapeAndGPmap

    void PrintLandscapeAndGPmapHeader(std::ostream& outfile){
        if(headerPrinted_==true) return;
        outfile<<"wAABB\twAABb\twAAbb\twAaBB\twAaBb\twAabb\twaaBB\twaaBb\twaabb";
        outfile<<"\tpAABB\tpAABb\tpAAbb\tpAaBB\tpAaBb\tpAabb\tpaaBB\tpaaBb\tpaabb"<<std::endl;
        headerPrinted_=true;
        }//PrintLandscapeAndGPmapHeader
    
    


    void ReportToScreen(void){//call coutLock before entering this function so it's not double-called
//        coutLock.lock();
        std::cout<<"  "<<solution_.genotype()<<"\tcode="<<solution_.hetcode()<<"\tphat=";
        if(pNeutral_){
                std::cout<<"neutral";}
            else{
                std::cout<<ROUND(phat_,decimalDigitsToRound);}
        std::cout<<"\tqhat=";
        if(qNeutral_){
                std::cout<<"neutral";}
            else{
                std::cout<<ROUND(qhat_,decimalDigitsToRound);}
        std::cout<<"\tphBar="<<ROUND(phBar_,decimalDigitsToRound)<<std::endl;
        std::cout<<"  solution mismatch pattern="<<solution_.mismatchPattern()<<"\tcode="<<solution_.mismatchHetcode()<<std::endl;
        std::cout<<"   from reference genotypes:"<<std::endl;
        for(int i=0;i<referenceGenotypes_.size();++i){
            std::cout<<"    "<<referenceGenotypes_[i].genotype()<<std::endl;}
        std::cout<<std::endl;
//        coutLock.unlock();
        }//reportToScreen

    
    
    void Reset(void){
        solution_.Clear();
        wBarMax_=phBar_=phat_=qhat_= -one;//error code
        pNeutral_=qNeutral_=false;
        TFdosageAlleles_.clear();TFprodAlleles_.clear();cisAlleles_.clear();
        fitnessLandscapesAndGPmaps_.clear();
        referenceGenotypes_.clear();
        }//Reset


     void SetFitnessLandscape(long double wBarMax,long double pBarAtMax,
                  long double pAABB,long double pAABb,long double pAAbb,
                  long double pAaBB,long double pAaBb,long double pAabb,
                  long double paaBB,long double paaBb,long double paabb,
                  long double wAABB,long double wAABb,long double wAAbb,
                  long double wAaBB,long double wAaBb,long double wAabb,
                 long double waaBB,long double waaBb,long double waabb){
        FitnessLandscapeParameters flp(wBarMax,pBarAtMax,pAABB,pAABb,pAAbb,pAaBB,pAaBb,pAabb,paaBB,paaBb,paabb,
                                        wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb);
         bool found=false;
         for(int f=0;f<fitnessLandscapesAndGPmaps_.size();++f){
             if(flp==fitnessLandscapesAndGPmaps_[f]){
                 found=true;}
             }//f
         if(!found){
             fitnessLandscapesAndGPmaps_.push_back(flp);}
        }//SetFitnessLandscape

     
    

    void SetSolution(SimplestRegPathIndividual& solutionIndiv,
                     long double wBarMax, long double phBar, long double phat, long double qhat,
                     bool pNeutral, bool qNeutral, std::string& hetcode, int bitstringLen){
        solution_.CollectData(solutionIndiv,bitstringLen);
        TFdosageAlleles_.clear();TFprodAlleles_.clear();cisAlleles_.clear();
        TFdosageAlleles_.push_back(solutionIndiv.TFdosage_[0]);
        TFprodAlleles_.push_back(solutionIndiv.TFproduct_[0]);
        if(solutionIndiv.TFdosage_[0]!=solutionIndiv.TFdosage_[1] || solutionIndiv.TFproduct_[0]!=solutionIndiv.TFproduct_[1]){
            TFdosageAlleles_.push_back(solutionIndiv.TFdosage_[1]);
            TFprodAlleles_.push_back(solutionIndiv.TFproduct_[1]);}
        cisAlleles_.push_back(solutionIndiv.cis_[0]);
       if(solutionIndiv.cis_[0]!=solutionIndiv.cis_[1]){
           cisAlleles_.push_back(solutionIndiv.cis_[1]);}
        wBarMax_=wBarMax;
        phBar_=phBar;
        phat_=phat;
        qhat_=qhat;
        pNeutral_=pNeutral;
        qNeutral_=qNeutral;
        trueHetcode_=solution_.trueHetcode_=hetcode;//overrides any neutral cases
        mismatchPattern_=solution_.mismatchPattern_=solutionIndiv.mismatchGenotype();
        mismatchHetcode_=solution_.mismatchHetcode_=solutionIndiv.mismatchHetType();//doesn't override mismatch  cases - should it?

        }//SetSolution
    
    
    friend int operator==(FitnessMaximumBitstringSolution& fd1, FitnessMaximumBitstringSolution& fd2);
    friend int operator!=(FitnessMaximumBitstringSolution& fd1, FitnessMaximumBitstringSolution& fd2);
    };//class fitnessMaximumBitstringSolution

const long double FitnessMaximumBitstringSolution::tol_=0.000001;

int operator==(FitnessMaximumBitstringSolution& fd1, FitnessMaximumBitstringSolution& fd2){
    // !! --> note that these can have different lists of fitness landscapes & reference genotypes;
    // it only compares the solution data !!
    if(fd1.solution_==fd2.solution_
       && ABS(fd1.wBarMax_-fd2.wBarMax_)<=fd1.tol_ && ABS(fd1.phBar_-fd2.phBar_)<=fd1.tol_
       && fd1.TFdosageAlleles_==fd2.TFdosageAlleles_ && fd1.TFprodAlleles_==fd2.TFprodAlleles_
       && fd1.cisAlleles_==fd2.cisAlleles_
       && ABS(fd1.phat_-fd2.phat_)<=fd1.tol_ && ABS(fd1.qhat_-fd2.qhat_)<=fd1.tol_
       && fd1.pNeutral_==fd2.pNeutral_ && fd1.qNeutral_==fd2.qNeutral_){
        return true;}
    return false;}

int operator!=(FitnessMaximumBitstringSolution& fd1, FitnessMaximumBitstringSolution& fd2){
    return !(fd1==fd2);}






class FitnessMaximaBitstringSolutions{//plural
public:
    std::vector<std::vector<FitnessMaximumBitstringSolution> > solutionsByHetCode_;
    long double wBarMax_;
    long double Popt_;
    long double omega_;//fitness variance
    bool headerPrinted_;
    static const long double tol_;
public:
    FitnessMaximaBitstringSolutions(void){
        wBarMax_=Popt_=omega_= (double)-one;//error code
        headerPrinted_=false;
        for(int h=0;h<26;++h){
            std::vector<FitnessMaximumBitstringSolution> hetcategory;
            solutionsByHetCode_.push_back(hetcategory);}
    }
    
    FitnessMaximaBitstringSolutions(const FitnessMaximaBitstringSolutions& fm){
        wBarMax_=fm.wBarMax_;
        Popt_=fm.Popt_;
        omega_=fm.omega_;
        solutionsByHetCode_=fm.solutionsByHetCode_;
        headerPrinted_=fm.headerPrinted_;
        }
    
    ~FitnessMaximaBitstringSolutions(void){
        solutionsByHetCode_.clear();}
    
    
    FitnessMaximaBitstringSolutions& operator=(FitnessMaximaBitstringSolutions& fm){
        wBarMax_=fm.wBarMax_;
        Popt_=fm.Popt_;
        omega_=fm.omega_;
        solutionsByHetCode_=fm.solutionsByHetCode_;
        headerPrinted_=fm.headerPrinted_;
        return *this;
        }

    void AddRefGtypeToSolution(int hetcode, int item, BitstringGenotypeData& refgtype){
        solutionsByHetCode_[hetcode][item].AppendReferenceGenotype(refgtype);
        }
    
    
    void AppendNewSolution(FitnessMaximumBitstringSolution& soln){
        int hetcode=soln.solutionHetcodeToInt();
        solutionsByHetCode_[hetcode].push_back(soln);
        }
    
    void ClearSolutions(void){
        for(int hc=0;hc<solutionsByHetCode_.size();++hc){
            solutionsByHetCode_[hc].clear();}
        }

    unsigned long items(void){
        unsigned long n=0;
        for(int h=0;h<26;++h){
            n += solutionsByHetCode_[h].size();}
        return n;}

    unsigned long items(int hetcode){
        return solutionsByHetCode_[hetcode].size();}

    
    void PrintSolutionTable(std::ostream& outfile, double Popt, double omega, int Ntf){
        PrintSolutionTableHeader(outfile);
        PrintSolutionTableData(outfile, Popt, omega, Ntf);
        }//PrintSolutionTable
    
    
    void PrintSolutionTableHeader(std::ostream& outfile){
        if(headerPrinted_==true) return;
        outfile<<"genotypes that maximize fitness:"<<std::endl;
        outfile<<"combinations with the highest fitness:"<<std::endl;
        outfile<<"Popt\tomega\tNtf\twBar\thet type\thet code\tphat\tqhat\tmean phenotype\tsolution genotype\t";
        outfile<<"mismatch genotype\tties\t1st ref g'type";
        outfile<<"\twAABB\twAABb\twAAbb\twAaBB\twAaBb\twAabb\twaaBB\twaaBb\twaabb";
        outfile<<"\tpAABB\tpAABb\tpAAbb\tpAaBB\tpAaBb\tpAabb\tpaaBB\tpaaBb\tpaabb"<<std::endl;
        headerPrinted_=true;
        }//PrintSolutionTableHeader
    
    
    void PrintSolutionTableData(std::ostream& outfile, double Popt, double omega, int Ntf){
        std::string tab("\t");
        for(int h=0;h<solutionsByHetCode_.size();++h){
            std::string c=intToHetcode(h);
            for(int i=0;i<solutionsByHetCode_[h].size();++i){
                FitnessMaximumBitstringSolution& s=solutionsByHetCode_[h][i];
                for(int j=0;j<s.fitnessLandscapesAndGPmaps_.size();++j){
                    FitnessLandscapeParameters& flp = s.fitnessLandscapesAndGPmaps_[j];
                    outfile<<Popt<<tab<<omega<<tab<<Ntf<<tab<<ROUND(wBarMax_,decimalDigitsToRound)<<tab<<h<<tab<<c<<tab;
                    if(!(s.pNeutral())){outfile<<ROUND(s.phat_,decimalDigitsToRound);}
                    outfile<<tab;
                    if(!(s.qNeutral())){outfile<<ROUND(s.qhat_,decimalDigitsToRound);}
                    outfile<<tab;
                    outfile<<ROUND(s.phBar_,decimalDigitsToRound)<<tab<<s.solutionGenotype();// <<tab<<s.ABgenotype();
                    outfile<<tab<<s.solutionMismatchPattern()<<tab<<s.referenceGenotypes_.size()<<tab;
                    outfile<<s.referenceGenotypes_[0].genotype()<<tab;
                    if(flp.wAABB_!=-one){outfile<<ROUND(flp.wAABB_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.wAABb_!=-one){outfile<<ROUND(flp.wAABb_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.wAAbb_!=-one){outfile<<ROUND(flp.wAAbb_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.wAaBB_!=-one){outfile<<ROUND(flp.wAaBB_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.wAaBb_!=-one){outfile<<ROUND(flp.wAaBb_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.wAabb_!=-one){outfile<<ROUND(flp.wAabb_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.waaBB_!=-one){outfile<<ROUND(flp.waaBB_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.waaBb_!=-one){outfile<<ROUND(flp.waaBb_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.waabb_!=-one){outfile<<ROUND(flp.waabb_,decimalDigitsToRound);}outfile<<tab;

                    if(flp.pAABB_!=-one){outfile<<ROUND(flp.pAABB_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.pAABb_!=-one){outfile<<ROUND(flp.pAABb_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.pAAbb_!=-one){outfile<<ROUND(flp.pAAbb_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.pAaBB_!=-one){outfile<<ROUND(flp.pAaBB_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.pAaBb_!=-one){outfile<<ROUND(flp.pAaBb_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.pAabb_!=-one){outfile<<ROUND(flp.pAabb_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.paaBB_!=-one){outfile<<ROUND(flp.paaBB_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.paaBb_!=-one){outfile<<ROUND(flp.paaBb_,decimalDigitsToRound);}outfile<<tab;
                    if(flp.paabb_!=-one){outfile<<ROUND(flp.paabb_,decimalDigitsToRound);}outfile<<std::endl;
                    }//j
                }//i
            }//h
        }//PrintSolutionTableData

    void ReplaceWbarMax(long double newWbarMax){
        solutionsByHetCode_.clear();
        wBarMax_=newWbarMax;//error code
        }
    
    void ReportToScreen(bool& newBest){//call coutLock before entering this function, so it's not double-called
//        coutLock.lock();
        std::cout<<"wBarMax="<<wBarMax_;
        if(newBest){std::cout<<"\t\tNEW BEST";newBest=false;}
        std::cout<<std::endl;
        for(int i=0;i<solutionsByHetCode_.size();++i){
            for(int j=0;j<solutionsByHetCode_[i].size();++j){
                solutionsByHetCode_[i][j].ReportToScreen();}}//i,j
        std::cout<<std::endl;
//        coutLock.unlock();
        }//ReportToScreen()
    
    
    void Reset(void){
        solutionsByHetCode_.clear();
        wBarMax_=Popt_=omega_= -one;//error code
        }
    
    void UpdateWithNewSolution(FitnessMaximumBitstringSolution& soln){
        //add solution if it has a novel bitstring genotype
        //if it's not, add it as a reference genotype to the relevant solution
        int hc = codeToInt(soln.trueHetcode_);
        if(soln.wBarMax_>wBarMax_){// || solutionsByHetCode_[hc].size()==0){
                ClearSolutions();
                solutionsByHetCode_[hc].push_back(soln);
                wBarMax_=soln.wBarMax_;
                }
            else if(soln.wBarMax_==wBarMax_){//
                bool found=false;
                for(unsigned long i=0;i<solutionsByHetCode_[hc].size();++i){
                    if(soln==solutionsByHetCode_[hc][i]){//checks for same bitstring only
                        found=true;
                        solutionsByHetCode_[hc][i].AppendReferenceGenotype(soln);
                        break;}
                    }//i
                if(!found){
                    AppendNewSolution(soln);}
                }
            else{}//don't add it
        }//UpdateWithNewSolution

    
    
    friend int operator==(FitnessMaximaBitstringSolutions& fm1, FitnessMaximaBitstringSolutions& fm2);
    friend int operator!=(FitnessMaximaBitstringSolutions& fm1, FitnessMaximaBitstringSolutions& fm2);

    };//class FitnessMaximaBitstringSolutions

const long double FitnessMaximaBitstringSolutions::tol_ = 0.000001;

int operator==(FitnessMaximaBitstringSolutions& fm1, FitnessMaximaBitstringSolutions& fm2){
    if(ABS(fm1.wBarMax_-fm2.wBarMax_)>fm1.tol_ || ABS(fm1.Popt_-fm2.Popt_)>fm1.tol_ || ABS(fm1.omega_-fm2.omega_)>fm1.tol_){
        return false;}
    for(unsigned long hc=0;hc<fm1.solutionsByHetCode_.size();++hc){
        if(fm1.solutionsByHetCode_[hc].size()!=fm2.solutionsByHetCode_[hc].size()){return false;}
        for(unsigned long i=0;i<fm1.solutionsByHetCode_[hc].size();++i){
            if(fm1.solutionsByHetCode_[hc][i]!=fm2.solutionsByHetCode_[hc][i]){ return false;}
            }}
    return true;
    }

int operator!=(FitnessMaximaBitstringSolutions& fm1, FitnessMaximaBitstringSolutions& fm2){
    return !(fm1==fm2);}



class FitnessMaximumSolutionSet{
    public:
    long double Popt_, omega_;
    long double Ntf_;
    long double wBarMax_;
    long double meanPhenotype_;
    long double p_,q_;//most-common allele frequencies
    bool pNeutral_,qNeutral_;
    int trueHetCode_;
    std::string trueHetPattern_;//3-character code
    int mismatchHetCode_;
    std::string mismatchHetPattern_;//3-character code
    std::string mismatchPattern_;//Mathematica-formatted code, e.g. {{0,0},{{0,0},{0,0}}}
    long numDuplicates_;//#reference genotypes giving this solution
    std::string firstSolutionGtype_;
    std::string firstRefGtype_;
    int bitstringLen_;
    bool splitSinglePoptRun_;
    uint64_t startingTF0val_, endTF0val_;
    static const long double tol_;

    public:
    FitnessMaximumSolutionSet(void):pNeutral_(false),qNeutral_(false),trueHetCode_(-1),
                                    mismatchHetCode_(-1),numDuplicates_(0),bitstringLen_(0),
                                    splitSinglePoptRun_(false),startingTF0val_(0),endTF0val_(0){
        Popt_=omega_=Ntf_=wBarMax_=meanPhenotype_=p_=q_=double(-one);}
                                    
    
    FitnessMaximumSolutionSet(int bitstringLen):pNeutral_(false),qNeutral_(false),trueHetCode_(-1),
                                    mismatchHetCode_(-1),numDuplicates_(0),bitstringLen_(bitstringLen),
                                    splitSinglePoptRun_(false),startingTF0val_(0){
        Popt_=omega_=Ntf_=wBarMax_=meanPhenotype_=p_=q_=double(-one);
        endTF0val_ = (uint64_t) pow(2,bitstringLen_)-1;}
                                    
    
    FitnessMaximumSolutionSet(int bitstringLen, bool splitSinglePoptRun, uint64_t startingTF0val, uint64_t endTF0val):
                                    pNeutral_(false),qNeutral_(false),trueHetCode_(-1),
                                    mismatchHetCode_(-1),numDuplicates_(0),bitstringLen_(bitstringLen),
                                    splitSinglePoptRun_(splitSinglePoptRun),startingTF0val_(startingTF0val),
                                    endTF0val_(endTF0val){
        Popt_=omega_=Ntf_=wBarMax_=meanPhenotype_=p_=q_=double(-one);}
                                    
    
    FitnessMaximumSolutionSet(long double Popt, long double omega, long double Ntf, int bitstringLen, long double wBarMax, long double meanPhenotype,
                                long double p, long double q, bool pNeutral, bool qNeutral,
                                int trueHetCode, std::string& trueHetPattern,
                                int mismatchHetCode, std::string& mismatchHetPattern,
                                std::string& mismatchPattern):numDuplicates_(1),
                                splitSinglePoptRun_(false),startingTF0val_(0){
        Popt_=Popt;
        omega_=omega;
        Ntf_=Ntf;
        wBarMax_=wBarMax;
        meanPhenotype_=meanPhenotype;
        p_=p;
        q_=q;
        pNeutral_=pNeutral;
        qNeutral_=qNeutral;
        trueHetCode_=trueHetCode;
        trueHetPattern_=trueHetPattern;
        mismatchHetCode_=mismatchHetCode;
        mismatchHetPattern_=mismatchHetPattern;
        mismatchPattern_=mismatchPattern;
        bitstringLen_=bitstringLen;
        endTF0val_ = (uint64_t) pow(2,bitstringLen_)-1;
        }

    FitnessMaximumSolutionSet(long double Popt, long double omega, long double Ntf, int bitstringLen, long double wBarMax, long double meanPhenotype,
                                long double p, long double q, bool pNeutral, bool qNeutral,
                                int trueHetCode, std::string& trueHetPattern,
                                int mismatchHetCode, std::string& mismatchHetPattern,
                                std::string& mismatchPattern,std::string& firstSolutionGtype,
                                std::string& firstRefGtype):numDuplicates_(1),splitSinglePoptRun_(false),startingTF0val_(0){
        Popt_=Popt;
        omega_=omega;
        Ntf_=Ntf;
        wBarMax_=wBarMax;
        meanPhenotype_=meanPhenotype;
        p_=p;
        q_=q;
        pNeutral_=pNeutral;
        qNeutral_=qNeutral;
        trueHetCode_=trueHetCode;
        trueHetPattern_=trueHetPattern;
        mismatchHetCode_=mismatchHetCode;
        mismatchHetPattern_=mismatchHetPattern;
        mismatchPattern_=mismatchPattern;
        firstSolutionGtype_=firstSolutionGtype;
        firstRefGtype_=firstRefGtype;
        bitstringLen_=bitstringLen;
        endTF0val_ = (uint64_t) pow(2,bitstringLen_)-1;
        }


    FitnessMaximumSolutionSet(long double Popt, long double omega, long double Ntf, int bitstringLen, long double wBarMax, long double meanPhenotype,
                                long double p, long double q, bool pNeutral, bool qNeutral,
                                int trueHetCode, std::string& trueHetPattern,
                                int mismatchHetCode, std::string& mismatchHetPattern,
                                std::string& mismatchPattern, bool splitSinglePoptRun,
                                uint64_t startingTF0val, uint64_t endTF0val):numDuplicates_(1),
                                splitSinglePoptRun_(splitSinglePoptRun),startingTF0val_(startingTF0val),
                                endTF0val_(endTF0val){
        Popt_=Popt;
        omega_=omega;
        Ntf_=Ntf;
        wBarMax_=wBarMax;
        meanPhenotype_=meanPhenotype;
        p_=p;
        q_=q;
        pNeutral_=pNeutral;
        qNeutral_=qNeutral;
        trueHetCode_=trueHetCode;
        trueHetPattern_=trueHetPattern;
        mismatchHetCode_=mismatchHetCode;
        mismatchHetPattern_=mismatchHetPattern;
        mismatchPattern_=mismatchPattern;
        bitstringLen_=bitstringLen;
        }

    FitnessMaximumSolutionSet(long double Popt, long double omega, long double Ntf, int bitstringLen, long double wBarMax, long double meanPhenotype,
                                long double p, long double q, bool pNeutral, bool qNeutral,
                                int trueHetCode, std::string& trueHetPattern,
                                int mismatchHetCode, std::string& mismatchHetPattern,
                                std::string& mismatchPattern,std::string& firstSolutionGtype,
                                std::string& firstRefGtype, bool splitSinglePoptRun,
                                uint64_t startingTF0val, uint64_t endTF0val):numDuplicates_(1),
                                splitSinglePoptRun_(splitSinglePoptRun),startingTF0val_(startingTF0val),
                                endTF0val_(endTF0val){
        Popt_=Popt;
        omega_=omega;
        Ntf_=Ntf;
        wBarMax_=wBarMax;
        meanPhenotype_=meanPhenotype;
        p_=p;
        q_=q;
        pNeutral_=pNeutral;
        qNeutral_=qNeutral;
        trueHetCode_=trueHetCode;
        trueHetPattern_=trueHetPattern;
        mismatchHetCode_=mismatchHetCode;
        mismatchHetPattern_=mismatchHetPattern;
        mismatchPattern_=mismatchPattern;
        firstSolutionGtype_=firstSolutionGtype;
        firstRefGtype_=firstRefGtype;
        bitstringLen_=bitstringLen;
        }


    FitnessMaximumSolutionSet(const FitnessMaximumSolutionSet& fmss){
        Popt_=fmss.Popt_;
        omega_=fmss.omega_;
        Ntf_=fmss.Ntf_;
        wBarMax_=fmss.wBarMax_;
        meanPhenotype_=fmss.meanPhenotype_;
        p_=fmss.p_;
        q_=fmss.q_;
        pNeutral_=fmss.pNeutral_;
        qNeutral_=fmss.qNeutral_;
        trueHetCode_=fmss.trueHetCode_;
        trueHetPattern_=fmss.trueHetPattern_;
        mismatchHetCode_=fmss.mismatchHetCode_;
        mismatchHetPattern_=fmss.mismatchHetPattern_;
        mismatchPattern_=fmss.mismatchPattern_;
        numDuplicates_=fmss.numDuplicates_;
        firstSolutionGtype_=fmss.firstSolutionGtype_;
        firstRefGtype_=fmss.firstRefGtype_;
        bitstringLen_=fmss.bitstringLen_;
        splitSinglePoptRun_=fmss.splitSinglePoptRun_;
        startingTF0val_=fmss.startingTF0val_;
        endTF0val_=fmss.endTF0val_;
        }
    
    ~FitnessMaximumSolutionSet(void){}
    
    FitnessMaximumSolutionSet& operator=(const FitnessMaximumSolutionSet& fmss){
        Popt_=fmss.Popt_;
        omega_=fmss.omega_;
        Ntf_=fmss.Ntf_;
        wBarMax_=fmss.wBarMax_;
        meanPhenotype_=fmss.meanPhenotype_;
        p_=fmss.p_;
        q_=fmss.q_;
        pNeutral_=fmss.pNeutral_;
        qNeutral_=fmss.qNeutral_;
        trueHetCode_=fmss.trueHetCode_;
        trueHetPattern_=fmss.trueHetPattern_;
        mismatchHetCode_=fmss.mismatchHetCode_;
        mismatchHetPattern_=fmss.mismatchHetPattern_;
        mismatchPattern_=fmss.mismatchPattern_;
        numDuplicates_=fmss.numDuplicates_;
        firstSolutionGtype_=fmss.firstSolutionGtype_;
        firstRefGtype_=fmss.firstRefGtype_;
        bitstringLen_=fmss.bitstringLen_;
        splitSinglePoptRun_=fmss.splitSinglePoptRun_;
        startingTF0val_=fmss.startingTF0val_;
        endTF0val_=fmss.endTF0val_;
        return *this;
        }

    public:
    long double Popt(void){return Popt_;}
    long double omega(void){return omega_;}
    long double Ntf(void){return Ntf_;}
    long double wBarMax(void){return wBarMax_;}
    long double meanPhenotype(void){return meanPhenotype_;}
    long double p(void){return p_;}
    long double q(void){return q_;}
    bool pNeutral(void){return pNeutral_;}
    bool qNeutral(void){return qNeutral_;}
    int trueHetCode(void){return trueHetCode_;}
    std::string& trueHetPattern(void){return trueHetPattern_;}
    int mismatchHetCode(void){return mismatchHetCode_;}
    std::string& mismatchHetPattern(void){return mismatchHetPattern_;}
    std::string& mismatchPattern(void){return mismatchPattern_;}
    long numDuplicates(void){return numDuplicates_;}
    std::string& firstSolutionGtype(void){return firstSolutionGtype_;}
    std::string& firstRefGtype(void){return firstRefGtype_;}

    void Increment(void){
        numDuplicates_++;}
    
    void PrintDataLine(std::ostream& outfile){
        std::string tab("\t");
        int Ntfint=(int)ROUND(Ntf_,0);
        outfile<<Popt_<<tab<<omega_<<tab<<bitstringLen_<<tab<<Ntfint<<tab<<wBarMax_<<tab<<meanPhenotype_<<tab;
        outfile<<trueHetPattern_<<tab<<trueHetCode_<<tab;
        if(!(pNeutral_)){outfile<<p_;} outfile<<tab;
        if(!(qNeutral_)){outfile<<q_;} outfile<<tab;
        outfile<<mismatchHetPattern_<<tab<<mismatchHetCode_<<tab<<mismatchPattern_<<tab;
        outfile<<firstSolutionGtype_<<tab<<firstRefGtype_<<tab<<numDuplicates_;
        if(splitSinglePoptRun_){
            outfile<<tab<<startingTF0val_<<tab<<endTF0val_;}
        outfile<<std::endl;
        }//PrintDataLine

    void PrintHeaderLine(std::ostream& outfile){
        outfile<<"Popt\tomega\tNtf\tbitstring len\twBarMax\tmean p'type\thet pattern\thet code\tp\tq";
        outfile<<"\tmismatch het pattern\tmismatch het code\tmismatch pattern";
        outfile<<"\t1st solution g'type\t1st ref g'type\titems";
        if(splitSinglePoptRun_){
            outfile<<"\tstarting TF0\tend TF0";}
        outfile<<std::endl;
        }//PrintHeaderLine

    void Reset(void){//Popt_, omega_ & bitstringLen_ unchanged
        wBarMax_=meanPhenotype_=p_=q_= (double)-one;
        pNeutral_=qNeutral_=false;
        trueHetCode_=mismatchHetCode_= -1;
        trueHetPattern_=mismatchHetPattern_=mismatchPattern_="";
        numDuplicates_=0;}
    
    void Set(long double Popt, long double omega, long double Ntf, int bitstringLen, long double wBarMax, long double meanPhenotype,
            long double p, long double q, bool pNeutral, bool qNeutral,
            int trueHetCode, std::string& trueHetPattern,
            int mismatchHetCode, std::string& mismatchHetPattern,
            std::string& mismatchPattern,std::string& firstSolutionGtype,
            std::string& firstRefGtype){
        Popt_=Popt;
        omega_=omega;
        Ntf_=Ntf;
        wBarMax_=wBarMax;
        meanPhenotype_=meanPhenotype;
        p_=p;
        q_=q;
        pNeutral_=pNeutral;
        qNeutral_=qNeutral;
        trueHetCode_=trueHetCode;
        trueHetPattern_=trueHetPattern;
        mismatchHetCode_=mismatchHetCode;
        mismatchHetPattern_=mismatchHetPattern;
        mismatchPattern_=mismatchPattern;
        firstSolutionGtype_=firstSolutionGtype;
        firstRefGtype_=firstRefGtype;
        numDuplicates_=1;
        bitstringLen_=bitstringLen;
        }

    void Set(long double Popt, long double omega, long double Ntf, int bitstringLen, long double wBarMax, long double meanPhenotype,
            long double p, long double q, bool pNeutral, bool qNeutral,
            int trueHetCode, std::string& trueHetPattern,
            int mismatchHetCode, std::string& mismatchHetPattern,
            std::string& mismatchPattern,std::string& firstSolutionGtype,
            std::string& firstRefGtype, bool splitSinglePoptRun,
            uint64_t startingTF0val, uint64_t endTF0val){
        Set(Popt,omega,Ntf,bitstringLen,wBarMax,meanPhenotype,p,q,pNeutral,qNeutral,trueHetCode,trueHetPattern,
                    mismatchHetCode,mismatchHetPattern,mismatchPattern,firstSolutionGtype,firstRefGtype);
        splitSinglePoptRun_=splitSinglePoptRun;
        startingTF0val_=startingTF0val;
        endTF0val_=endTF0val;
        }


    friend int operator==(const FitnessMaximumSolutionSet& fmss1, const FitnessMaximumSolutionSet& fmss2);
    friend int operator!=(const FitnessMaximumSolutionSet& fmss1, const FitnessMaximumSolutionSet& fmss2);
    };//FitnessMaximumSolutionSet

const long double FitnessMaximumSolutionSet::tol_=0.000001;

int operator==(const FitnessMaximumSolutionSet& fmss1, const FitnessMaximumSolutionSet& fmss2){
    //doesn't compare bitstring length, first solution or reference gtypes
    if(fmss1.Popt_==fmss2.Popt_ && fmss1.omega_==fmss2.omega_ && fmss1.Ntf_==fmss2.Ntf_
       && ABS(fmss1.wBarMax_-fmss2.wBarMax_)<=fmss1.tol_
       && ABS(fmss1.meanPhenotype_-fmss2.meanPhenotype_)<=fmss1.tol_
       && ABS(fmss1.p_-fmss2.p_)<=fmss1.tol_ && ABS(fmss1.q_-fmss2.q_)<=fmss1.tol_
       && fmss1.pNeutral_==fmss2.pNeutral_ && fmss1.qNeutral_==fmss2.qNeutral_
       && fmss1.trueHetCode_==fmss2.trueHetCode_ && fmss1.trueHetPattern_==fmss2.trueHetPattern_
       && fmss1.mismatchHetCode_==fmss2.mismatchHetCode_ && fmss1.mismatchHetPattern_==fmss2.mismatchHetPattern_
       && fmss1.mismatchPattern_==fmss2.mismatchPattern_ && fmss1.splitSinglePoptRun_==fmss2.splitSinglePoptRun_
       && fmss1.startingTF0val_==fmss2.startingTF0val_ && fmss1.endTF0val_==fmss2.endTF0val_){return true;}
    return false;
    }

int operator!=(const FitnessMaximumSolutionSet& fmss1, const FitnessMaximumSolutionSet& fmss2){
    return !(fmss1==fmss2);}




class FitnessMaximaSolutionSets{
    public:
    std::vector<long double> wBarMaxPerPopt_;
    std::vector<long double> PoptValuesStored_;
    std::vector<std::vector<FitnessMaximumSolutionSet> > uniqueSolutionsByPopt_;//[Popt][solution]
    bool splitSinglePoptRun_;
    uint64_t startingTF0val_, endTF0val_;
    
    public:
    FitnessMaximaSolutionSets(void):splitSinglePoptRun_(false),startingTF0val_(0),endTF0val_(0){}

    FitnessMaximaSolutionSets(bool splitSinglePoptRun, uint64_t startingTF0val, uint64_t endTF0val):
            splitSinglePoptRun_(splitSinglePoptRun),startingTF0val_(startingTF0val),endTF0val_(endTF0val){}

    FitnessMaximaSolutionSets(const FitnessMaximaSolutionSets& fms){
        wBarMaxPerPopt_=fms.wBarMaxPerPopt_;
        PoptValuesStored_=fms.PoptValuesStored_;
        uniqueSolutionsByPopt_=fms.uniqueSolutionsByPopt_;
        splitSinglePoptRun_=fms.splitSinglePoptRun_;
        startingTF0val_=fms.startingTF0val_;
        endTF0val_=fms.endTF0val_;}
        
    ~FitnessMaximaSolutionSets(void){}
    
    FitnessMaximaSolutionSets& operator=(FitnessMaximaSolutionSets& fms){
        wBarMaxPerPopt_=fms.wBarMaxPerPopt_;
        PoptValuesStored_=fms.PoptValuesStored_;
        uniqueSolutionsByPopt_=fms.uniqueSolutionsByPopt_;
        splitSinglePoptRun_=fms.splitSinglePoptRun_;
        startingTF0val_=fms.startingTF0val_;
        endTF0val_=fms.endTF0val_;
        return *this;}
    
    public:
    
    void AddSolution(FitnessMaximumSolutionSet& fmss){
        if(!(includesSolutionsForPopt(fmss.Popt()))){//make a new Popt and store this solution
            wBarMaxPerPopt_.push_back(fmss.wBarMax());
            PoptValuesStored_.push_back(fmss.Popt());
            std::vector<FitnessMaximumSolutionSet> v;
            v.push_back(fmss);
            uniqueSolutionsByPopt_.push_back(v);
            return;}
        int p=indexForPopt(fmss.Popt());
        long double wMaxNew = fmss.wBarMax(), wMaxOld = uniqueSolutionsByPopt_[p][0].wBarMax();
        if(wMaxNew<wMaxOld) return;//ignore
        if(wMaxNew>wMaxOld){//clear old values
            wBarMaxPerPopt_[p]=wMaxNew;
            uniqueSolutionsByPopt_[p].clear();
            uniqueSolutionsByPopt_[p].push_back(fmss);
            return;}
        bool found=false;
        for(int i=0;i<uniqueSolutionsByPopt_[p].size();++i){
            if(fmss==uniqueSolutionsByPopt_[p][i]){//increment if it's already on the books
                uniqueSolutionsByPopt_[p][i].Increment();
                found=true;}
            if(found)break;
            }//i
        if(!found){//it's unique, so add it
            uniqueSolutionsByPopt_[p].push_back(fmss);}
        }//AddSolution


    void ConcatenateSolutions(FitnessMaximaSolutionSets& fmsSets){
        for(int p=0;p<fmsSets.uniqueSolutionsByPopt_.size();++p){
            for(int fmss=0;fmss<fmsSets.uniqueSolutionsByPopt_[p].size();++fmss){
                AddSolution(fmsSets.uniqueSolutionsByPopt_[p][fmss]);
                }//fmss
            }//p
        }//ConcatenateSolutions


    bool includesSolutionsForPopt(long double Popt){
        for(unsigned long i=0;i<PoptValuesStored_.size();++i){
            if(PoptValuesStored_[i]==Popt){ return true;}}
        return false;}

    int indexForPopt(long double Popt){//returns out of range error if missing (-1)
        for(int i=0;i<PoptValuesStored_.size();++i){
            if(PoptValuesStored_[i]==Popt){ return i;}}
        return -1;}
    
    void PrintDataByPopt(std::ostream& outfile, long double Popt){
        if(!(includesSolutionsForPopt(Popt)))return;
        int p = indexForPopt(Popt);
        for(int i=0;i<uniqueSolutionsByPopt_[p].size();++i){
            uniqueSolutionsByPopt_[p][i].PrintDataLine(outfile);}
        }//PrintDataByPopt


    void PrintHeaderLine(std::ostream& outfile){
        outfile<<"Popt\tomega\tbitstring len\tNtfsat\twBarMax\tmean phenotype\thet pattern\thet code\tp\tq";
        outfile<<"\tmismatch het pattern\tmismatch het code\tmismatch pattern\t";
        outfile<<"1st solution g'type\t1st ref g'type\titems";
        if(splitSinglePoptRun_){
            outfile<<"\tstarting TF0\tend TF0";}
        outfile<<std::endl;
        }

    long double wBarMax(long double Popt){
        int p = indexForPopt(Popt);
        if(p==-1) {return -one;}
        return wBarMaxPerPopt_[p];}
    
    
    friend int operator==(FitnessMaximaSolutionSets& fms1, FitnessMaximaSolutionSets& fms2);
    friend int operator!=(FitnessMaximaSolutionSets& fms1, FitnessMaximaSolutionSets& fms2);
    };// class FitnessMaximaSolutionSets





int operator==(FitnessMaximaSolutionSets& fms1, FitnessMaximaSolutionSets& fms2){
    return (fms1.wBarMaxPerPopt_==fms2.wBarMaxPerPopt_ && fms1.PoptValuesStored_==fms2.PoptValuesStored_
            && fms1.uniqueSolutionsByPopt_==fms2.uniqueSolutionsByPopt_);}

int operator!=(FitnessMaximaSolutionSets& fms1, FitnessMaximaSolutionSets& fms2){
    return !(fms1==fms2);}









long double meanPhenotype(long double ph_AABB,long double ph_AABb,long double ph_AAbb,
                 long double ph_AaBB,long double ph_AaBb,long double ph_Aabb,
                 long double ph_aaBB,long double ph_aaBb,long double ph_aabb,
                 long double p, long double q){
    //p is the freq of TF allele A; q is the frequency of cis allele B
    long double mean = ph_AABB*p*p*q*q + two*ph_AABb*p*p*q*(one - q) + ph_AAbb*p*p*(one - q)*(one - q)
              + two*ph_AaBB*p*(one - p)*q*q + four*ph_AaBb*p*(one - p)*q*(one - q)
              + two*ph_Aabb*p*(one - p)*(one - q)*(one - q)
              + ph_aaBB*(one - p)*(one - p)*q*q + two*ph_aaBb*(one - p)*(one - p)*q*(one - q)
              + ph_aabb*(one - p)*(one - p)*(one - q)*(one - q);

    return mean;
    }//meanPhenotype


long double meanPhenotype(long double ph_AA,long double ph_Aa,long double ph_aa,long double p){
    //p is the frequency of allele A
    long double mean = ph_AA*p*p + two*ph_Aa*p*(one - p) + ph_aa*(one - p)*(one - p);
    return mean;
    }


void collectMeanPhenotypesPQ(long double ph_AABB,long double ph_AABb,long double ph_AAbb,
                            long double ph_AaBB,long double ph_AaBb,long double ph_Aabb,
                            long double ph_aaBB,long double ph_aaBb,long double ph_aabb,
                            std::vector<long double>& phat, std::vector<long double>& qhat,
                            std::vector<long double>& meanPhenotypes, int& numMaxima){
    for(int i=0;i<meanPhenotypes.size();++i){meanPhenotypes[i]=-one;}
    for(int m=0;m<numMaxima;++m){
        meanPhenotypes[m]=meanPhenotype(ph_AABB,ph_AABb,ph_AAbb,ph_AaBB,ph_AaBb,ph_Aabb,ph_aaBB,ph_aaBb,ph_aabb,
                            MIN(MAX(zero,phat[m]),one), MIN(MAX(zero,qhat[m]),one));//accounts for error flag when p or q is neutral
        }//m
    }//collectMeanPhenotypesPQ


void collectMeanPhenotypesP(long double ph_AA,long double ph_Aa,long double ph_aa,
                            std::vector<long double>& phat, std::vector<long double>& meanPhenotypes, int& numMaxima){
    for(int m=0;m<numMaxima;++m){
        meanPhenotypes[m]=meanPhenotype(ph_AA,ph_Aa,ph_aa,phat[m]);
        }//m
    }//collectMeanPhenotypesP

void collectMeanPhenotypesQ(long double ph_BB,long double ph_Bb,long double ph_bb,
                            std::vector<long double>& qhat, std::vector<long double>& meanPhenotypes, int& numMaxima){
    for(int m=0;m<numMaxima;++m){
        meanPhenotypes[m]=meanPhenotype(ph_BB,ph_Bb,ph_bb,qhat[m]);
        }//m
    }//collectMeanPhenotypesQ




long double wBar(long double wAABB,long double wAABb,long double wAAbb,
				 long double wAaBB,long double wAaBb,long double wAabb,
				 long double waaBB,long double waaBb,long double waabb,
				 long double p, long double q){
	//p is the freq of TF allele A; q is the frequency of cis allele B
	long double wbar = wAABB*p*p*q*q + two*wAABb*p*p*q*(one - q) + wAAbb*p*p*(one - q)*(one - q)
  			+ two*wAaBB*p*(one - p)*q*q + four*wAaBb*p*(one - p)*q*(one - q)
  			+ two*wAabb*p*(one - p)*(one - q)*(one - q)
  			+ waaBB*(one - p)*(one - p)*q*q + two*waaBb*(one - p)*(one - p)*q*(one - q)
  			+ waabb*(one - p)*(one - p)*(one - q)*(one - q);

	return wbar;
	}

long double wBar(long double wAA,long double wAa,long double waa,long double p){
	//p is the frequency of allele A
	long double wbar = wAA*p*p + two*wAa*p*(one - p) + waa*(one - p)*(one - p);
	return wbar;
	}



long double DwBar_dp(long double wAABB,long double wAABb,long double wAAbb,
				 long double wAaBB,long double wAaBb,long double wAabb,
				 long double waaBB,long double waaBb,long double waabb,
				 long double p, long double q){
	//p is the frequency of TF allele A; q is the frequency of cis allele B
/*	//derivative solved using Mathematica
D[wbar[wAABB, wAABb, wAAbb, wAaBB, wAaBb, wAabb, waaBB, waaBb, waabb,
p, q], p] // Simplify
   	= 2 ((-1 + p) (-1 + q)^2 waabb + wAabb - 2 p wAabb + p wAAbb +
	2 q ((-1 + p) waaBb + (-1 + 2 p) wAabb + wAaBb - 2 p wAaBb -
	   p wAAbb + p wAABb) +
	q^2 (-2 (-1 + p) waaBb + (-1 + p) waaBB + wAabb - 2 p wAabb -
	   2 wAaBb + 4 p wAaBb + wAaBB - 2 p wAaBB + p wAAbb - 2 p wAABb +
	   p wAABB))
   */

	long double d =
   	two*((-one + p)*(-one + q)*(-one + q)*waabb + wAabb - two*p*wAabb + p*wAAbb +
	two*q*((-one + p)*waaBb + (-one + two*p)*wAabb + wAaBb - two*p*wAaBb -
	   p*wAAbb + p*wAABb) +
	q*q*(-two*(-one + p)*waaBb + (-one + p)*waaBB + wAabb - two*p*wAabb -
	   two*wAaBb + four*p*wAaBb + wAaBB - two*p*wAaBB + p*wAAbb - two*p*wAABb +
	   p*wAABB));

	return d;
	}//DwBar_dp


long double DwBar_dp(long double wAA, long double wAa, long double waa, long double p){
	//p is the frequency of TF allele A
/*	//derivative solved using Mathematica
	D[wAA p^2 + 2 wAa p (1 - p) + waa (1 - p)^2, p] // Simplify
		= 2 (wAa - waa + p (wAA - 2 wAa + waa))
*/
	long double d = two*(wAa - waa + p*(wAA - two*wAa + waa));
	return d;
	}//DwBar_dp


long double DwBar_dq(long double wAABB,long double wAABb,long double wAAbb,
				 long double wAaBB,long double wAaBb,long double wAabb,
				 long double waaBB,long double waaBb,long double waabb,
				 long double p, long double q){
	//p is the frequency of TF allele A; q is the frequency of cis allele B
/*
	//derivative solved using Mathematica
D[wAABB p^2 q^2 + 2 wAABb p^2 q (1 - q) + wAAbb p^2 (1 - q)^2 +
 2 wAaBB p (1 - p) q^2 + 4 wAaBb p (1 - p) q (1 - q) +
 2 wAabb  p (1 - p) (1 - q)^2 + waaBB (1 - p)^2 q^2 +
 2 waaBb (1 - p)^2*q*(1 - q) + waabb (1 - p)^2 (1 - q)^2 ,
q] // Simplify
  = 2 ((-1 + p)^2 (-1 + q) waabb - (-1 + p)^2 (-1 + q) waaBb - (-1 +
	 p)^2 q waaBb + (-1 + p)^2 q waaBB -
  2 (-1 + p) p (-1 + q) wAabb + 2 (-1 + p) p (-1 + q) wAaBb +
  2 (-1 + p) p q wAaBb - 2 (-1 + p) p q wAaBB + p^2 (-1 + q) wAAbb -
  p^2 (-1 + q) wAABb - p^2 q wAABb + p^2 q wAABB)
*/
	long double d =
		two*((-one + p)*(-one + p)*(-one + q)*waabb - (-one + p)*(-one + p)*(-one + q)*waaBb - (-one +
		   p)*(-one + p)*q*waaBb + (-one + p)*(-one + p)*q*waaBB -
		two*(-one + p)*p*(-one + q)*wAabb + two*(-one + p)*p*(-one + q)*wAaBb +
		two*(-one + p)*p*q*wAaBb - two*(-one + p)*p*q*wAaBB + p*p*(-one + q)*wAAbb -
		p*p*(-1 + q)*wAABb - p*p*q*wAABb + p*p*q*wAABB);
	 return d;
	 }//DwBar_dq


long double DwBar_dq(long double wBB,long double wBb,long double wbb,long double q){
	//q is the frequency of cis allele B
/*	derivative solved using Mathematica
	D[wBB q^2 + 2 wBb q (1 - q) + wbb (1 - q)^2, q] // Simplify
		= 2 (wBb - wbb + q (wBB - 2 wBb + wbb))
*/
	long double d = two*(wBb - wbb + q*(wBB - two*wBb + wbb));
	 return d;
	 }//DwBar_dq

//******************* Press et al. 1992 algorithms for function minimization/maximization ***********
//a Press et al. 1992 algorithm used in 1D function maximization, modified for maximizing wBar with p

inline void Shift(long double& a, long double& b, long double& c, long double& d){
	a=b;b=c;c=d;}


/* **************** modified From Numerical Recipes in C **********************/

void get_psum(long double *pABC, long double *qABC, long double& psum, long double& qsum){
	psum=qsum=zero;
	for (int i=0;i<3;i++){
		psum += pABC[i]; qsum += qABC[i];}
	}//get_psum


void midpoint(long double pA, long double qA, long double pB, long double qB,
					long double& pMidpoint, long double& qMidpoint){
	//distance from pt C to the midpoint between points A & B
	pMidpoint= (pA+pB)/two;
	qMidpoint= (qA+qB)/two;
	}//midpoint


void distanceToMidpoint(long double pC, long double qC, long double pA, long double qA, long double pB, long double qB,
					long double& pMidpointDist, long double& qMidpointDist){
	//distance from pt C to the midpoint between points A & B
	pMidpointDist= ((pA-pC)+(pB-pC))/two;
	qMidpointDist= ((qA-qC)+(qB-qC))/two;
	}//distanceToMidpoint


void SortByWbar(long double *wbarABC, long double *pABC, long double *qABC){
	if(wbarABC[0]>=wbarABC[1] && wbarABC[1]>=wbarABC[2]){//already sorted
		return;}
	//This sorts lists coordinates {p,q} at points A, B and C on the fitness landscape,
	//in the order of their mean fitnesses
	//orders can be:
	//a) high,med,low --> correct order
	//b) high,low,med
	//c) med,high,low
	//d) med,low,high
	//e) low,high,med
	//f) low,med,high
	if(wbarABC[1]>wbarABC[0]){
		//e) low,high,med --> b) high,low,med
		//f) low,med,high --> d) med,low,high
		//c) med,high,low --> a) high,med,low --> (c->a) done
		SWAP(wbarABC[1],wbarABC[0]);
		SWAP(pABC[1],pABC[0]);
		SWAP(qABC[1],qABC[0]);}
	if(wbarABC[2]>wbarABC[1]){
		//b) high,low,med --> a) high,med,low --> 	(b->a && e->b->a) done
		//d) med,low,high --> c) med,high,low -->	(d->c && f->d->c)
		SWAP(wbarABC[2],wbarABC[1]);
		SWAP(pABC[2],pABC[1]);
		SWAP(qABC[2],qABC[1]);}
	if(wbarABC[1]>wbarABC[0]){
		//c) med,high,low --> a) high,med,low --> (d->c->a && f->d->c->a) done
		SWAP(wbarABC[1],wbarABC[0]);
		SWAP(pABC[1],pABC[0]);
		SWAP(qABC[1],qABC[0]);}
	}//SortByWbar()




bool wBarMaximizedToWithinTolerance(long double *wbarABC, long double *pABC, long double *qABC, long double tol){
	//true if the fitnesses and points are sufficiently close together
	if(wbarABC[0]-wbarABC[2]<=tol){
		if(ABS(pABC[0]-pABC[1])<=tol && ABS(pABC[0]-pABC[2])<=tol && ABS(pABC[1]-pABC[2])<=tol){
			if(ABS(qABC[0]-qABC[1])<=tol && ABS(qABC[0]-qABC[2])<=tol && ABS(qABC[1]-qABC[2])<=tol){
				return true;}}}
	return false;
	}//wBarMaximizedToWithinTolerance






long double amoebaTry(long double w0000,long double w0001,long double w0011,
								 long double w0100,long double w0101,long double w0111,
								 long double w1100,long double w1101,long double w1111,
								 long double wbarAtC, long double& pAtC, long double& qAtC,
								 long double& pAtABmidpoint, long double& qAtABmidpoint,//RanNumGen ran,
								 long double factor){
	//modified from Press et al.'s amtry() -- follows Wikipedia's description of the algorithm (https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method)
	//Extrapolates a random distance, scaled by a factor fac, through the face of the simplex across from the minimum,
	//tries it, and replaces the minimum point if the new point is better.
	long double ptry,qtry;
	ptry=pAtABmidpoint + factor*(pAtABmidpoint-pAtC);
	qtry=qAtABmidpoint + factor*(qAtABmidpoint-qAtC);
	ptry=MIN(MAX(ptry,zero),one);
	qtry=MIN(MAX(qtry,zero),one);
	//Evaluate the function at the trial point.
	long double wbarTry=wBar(w0000,w0001,w0011,w0100,w0101,w0111,w1100,w1101,w1111,ptry,qtry);
	if (wbarTry > wbarAtC){// If it’s better than the lowest, then replace the highest.
		pAtC=ptry;
		qAtC=qtry;
		}
	return wbarTry;
	}//amoebaTry




int amoeba(long double w0000,long double w0001,long double w0011,
								 long double w0100,long double w0101,long double w0111,
								 long double w1100,long double w1101,long double w1111,
								 long double *wbarABC, long double *pABC, long double *qABC,
								 //points A, B & C must be sorted must be sorted by fitness from high-low
								 long double tolerance, int& nFuntionCalls){
/*
Multidimensional maximization of wBar() by the inverting the downhill simplex method of Nelder and Mead.
Basically, it takes points A, B and C and swaps them around until they converge on a maximum, within a
predetermined margin of error.  Point A is always set to the maximum as the points are shifted around.
wbarABC[0] is the maximum fitness, and the coordinates are {pABC[0],qABC[0]}.  It maximizes wBar in the
center of the {p,q} space, and the function that calls it checks that against the edges & corners to
get the final maximum.
*/
	int maxRepeatsWithoutImprovement=100;
	nFuntionCalls=0;
	int repeatsWithoutImprovement=0, iterations=0;
	bool maximized=false,betterPointFound=false;
	long double pAtABmidpoint,qAtABmidpoint;
	long double wbarlow;//,wbarhigh;
	long double wbarToTry1,wbarToTry2,wbarToTry3;
	long double plow1,qlow1,plow2,qlow2,plow3,qlow3;//pbest,plow,qbest,qlow
	long double factor1=one,factor2=two,factor3=half;
	while(!maximized){
		iterations++;betterPointFound=false;
		midpoint(pABC[0],qABC[0],pABC[1],qABC[1],pAtABmidpoint,qAtABmidpoint);
		plow1=plow2=plow3=pABC[2];
		qlow1=qlow2=qlow3=qABC[2];wbarlow=wbarABC[2];
		//try from C partway to the AB midpoint,
		//then from C away from the AB midpoint,
		//then from C to beyond the AB midpoint
	wbarToTry1=amoebaTry(w0000,w0001,w0011,w0100,w0101,w0111,w1100,w1101,w1111,wbarlow,plow1,qlow1,pAtABmidpoint,qAtABmidpoint,factor1);
		if(wbarToTry1<wbarABC[0] && wbarToTry1>wbarABC[1]){//if it's the 2nd best
				wbarABC[2]=wbarToTry1;pABC[2]=plow1;qABC[2]=qlow1;betterPointFound=true;//replace the last one with it and re-sort
				}
			else if(wbarToTry1>wbarABC[0]){//if it's the best, then go further using it as the starting point
				plow2=plow1;qlow2=qlow1;
				wbarToTry2=amoebaTry(w0000,w0001,w0011,w0100,w0101,w0111,w1100,w1101,w1111,wbarlow,plow2,qlow2,
						pAtABmidpoint,qAtABmidpoint,factor2);
				if(wbarToTry2>wbarToTry1){//if it's even better, use it instead
						wbarABC[2]=wbarToTry2;pABC[2]=plow2;qABC[2]=qlow2;
						betterPointFound=true;repeatsWithoutImprovement=0;}
					else{
						wbarABC[2]=wbarToTry1;pABC[2]=plow1;qABC[2]=qlow1;
						betterPointFound=true;repeatsWithoutImprovement=0;}
				}//wbarToTry1>wbarABC[0]
			else{// if(wbarToTry1<wbarABC[1]){//but maybe not wbarToTry1<wbarABC[2]
				wbarToTry3=amoebaTry(w0000,w0001,w0011,w0100,w0101,w0111,w1100,w1101,w1111,wbarlow,plow3,qlow3,
						pAtABmidpoint,qAtABmidpoint,factor3);
				if(wbarToTry3>wbarABC[2]){
					wbarABC[2]=wbarToTry3;pABC[2]=plow3;qABC[2]=qlow3;
					betterPointFound=true;repeatsWithoutImprovement=0;}
				}
		if(betterPointFound){
				SortByWbar(wbarABC,pABC,qABC);factor1=one;factor3=half;
				maximized=wBarMaximizedToWithinTolerance(wbarABC,pABC,qABC,tolerance);}
			else{
				repeatsWithoutImprovement++;
				factor1 *=(long double)0.9;factor3*=(long double)0.9;
				if(factor1<=one/ten){//sufficiently close I guess
					maximized=true;}
				}
		if(!maximized && repeatsWithoutImprovement>=maxRepeatsWithoutImprovement){
            coutLock.lock();
            std::cout<<"can\'t maximize fitness for case where";
			std::cout<<" w0000="<<w0000;
			std::cout<<"; w0001="<<w0001;
			std::cout<<"; w0011="<<w0011;
			std::cout<<"; w0100="<<w0100;
			std::cout<<"; w0100="<<w0101;
			std::cout<<"; w0100="<<w0111;
			std::cout<<"; w1100="<<w1100;
			std::cout<<"; w1100="<<w1101;
			std::cout<<"; w1100="<<w1111;
			std::cout<<std::endl;
			std::cout<<"Skipping this case"<<std::endl;
            coutLock.unlock();
			return 0;
			}//too many repeats without improvement
		}//while !maximized
	return 1;
	}//end function amoeba

/* *************************************************** */


void MaximizePopMeanFitnessP(long double wAA, long double wAa, long double waa,
							 long double& wBarMax, std::vector<long double>& phat,
                             std::vector<bool>& pNeutral, int& numMaxima){
	//at wBarMax, phat[] holds the equilibrium frequency/frequencies of allele A
		//there can be 2 if wAA==waa>wAa
	//in the notation below, P0 is allele frequency at point p=0 = 'a' fixed
	//						P1 is the frequency at point p=1 = 'A' fixed
	numMaxima=1;
    pNeutral[0]=false;
	if(wAA==waa){//easy one
		if(wAa>wAA){//heterozygote advantage
				long double pp;
				phat[0]=pp=half;
				wBarMax = wBar(wAA,wAa,waa,pp);}
			else if(wAA>wAa){//heterozygote disadvantage; symmetrical, unstable equilibrium
				numMaxima=2;
                phat[0]=one;phat[1]=zero;pNeutral[1]=false;wBarMax=wAA;}
			else{//completely neutral; set to 0.5 arbitrarily
                phat[0]=phat[1]=half; pNeutral[0]=true; wBarMax=wAA;}
		return;}//wAA==waa
	long double dwBarAtP0 = DwBar_dp(wAA,wAa,waa,zero);//a fixed
	long double dwBarAtP1 = DwBar_dp(wAA,wAa,waa,one);//A fixed
	if(dwBarAtP0<zero && dwBarAtP1>zero){//wAa is lowest; no polymorphism maintained
			if(wAA>=waa){
					wBarMax=wAA;phat[0]=one;}//if they're equal, this sets p=1 arbitrarily
				else{
					wBarMax=waa;phat[0]=zero;}
			}
		else if(dwBarAtP0>zero && dwBarAtP1>zero){//steadily climbing towards p=1
			wBarMax=wAA;phat[0]=one;}
		else if(dwBarAtP0<zero && dwBarAtP1<zero){//steadily climbing towards p=0
			wBarMax=waa;phat[0]=zero;}
		else{//highest in the middle
			if(ABS(wAA-two*wAa+waa)==zero){phat[0]=one;}
				else{phat[0]=ABS((waa-wAa)/(wAA-two*wAa+waa));}
			wBarMax=wBar(wAA,wAa,waa,phat[0]);
			}
	wBarMax=MIN(MAX(ROUND(wBarMax,decimalDigitsToRound),zero),one);//round & bound
	for(int m=0;m<numMaxima;++m){
		phat[m]=MIN(MAX(ROUND(phat[m],decimalDigitsToRound),zero),one);}
	}//MaximizePopMeanFitnessP


void MaximizePopMeanFitnessQ(long double wBB, long double wBb, long double wbb,
								 long double& wBarMax, std::vector<long double>& qhat,
                             std::vector<bool>& qNeutral, int& numMaxima){
	numMaxima=1;
    qNeutral[0]=false;
	if(wBB==wbb){//easy one
		if(wBb>wBB){//heterozygote advantage
				long double qq;
				qhat[0]=qq=half;
				wBarMax = wBB*qq*qq+wBb*two*qq*(one-qq)+wbb*(one-qq)*(one-qq);}
			else if(wBB>wBb){//heterozygote disadvantage; symmetrical, unstable equilibrium
				numMaxima=2;
                qNeutral[1]=false;
                qhat[0]=one;qhat[1]=zero;wBarMax=wBB;}
			else{//completely neutral; set to 0.5 arbitrarily
                qNeutral[0]=true; qhat[0]=half; wBarMax=wBb;}
		return;}
	long double dwBarAtQ0 = DwBar_dq(wBB,wBb,wbb,zero);
	long double dwBarAtQ1 = DwBar_dq(wBB,wBb,wbb,one);
	if(dwBarAtQ0<zero && dwBarAtQ1>zero){//no polymorphism maintained
			if(wBB>=wbb){
                    wBarMax=wBB;qhat[0]=one;}//if they're equal, this sets q=1 arbitrarily
				else{
                    wBarMax=wbb;qhat[0]=zero;}
			}
		else if(dwBarAtQ0>zero && dwBarAtQ1>zero){//steadily climbing towards q=1
			wBarMax=wBB;qhat[0]=one;}
		else if(dwBarAtQ0<zero && dwBarAtQ1<zero){//steadily climbing towards q=0
			wBarMax=wbb;qhat[0]=zero;}
		else{//highest in the middle
			if(ABS(wBB-two*wBb+wbb)==zero){qhat[0]=one;}
				else{qhat[0]=ABS((wbb-wBb)/(wBB-two*wBb+wbb));}
			wBarMax=wBar(wBB,wBb,wbb,qhat[0]);
			}
	wBarMax=MIN(MAX(ROUND(wBarMax,decimalDigitsToRound),zero),one);//round & bound
	for(int m=0;m<numMaxima;++m){
		qhat[m]=MIN(MAX(ROUND(qhat[m],decimalDigitsToRound),zero),one);}
	}//MaximizePopMeanFitnessQ






int MaximizePopMeanFitnessPandQv2(long double wAABB,long double wAABb,long double wAAbb,
								 long double wAaBB,long double wAaBb,long double wAabb,
								 long double waaBB,long double waaBb,long double waabb,long double& wBarMax,
                                  std::vector<long double>& phat, std::vector<long double>& qhat,
                                  std::vector<bool>& pNeutral, std::vector<bool>& qNeutral,
                                  int& numMaxima, std::string& focalGtypeStr){
	//find wBar at the 4 corners, then maximize wBar starting at p[A]=q[B]=0.5.  Finally, compare all 5 to get maximum
	//need to do this because fitness surface can be convex in the middle, but still highest in one of the corners
	//returns 0 if maximization fails; 1 otherwise
    //error codes:  phat=2 means p is neutral; qhat=2 means q is neutral

	//initialize

	//start at each corner
    for(int i=0;i<phat.size();++i){phat[i]=qhat[i]=-one;}
    for(int i=0;i<pNeutral.size();++i){pNeutral[i]=qNeutral[i]=false;}
    bool solutionFound=false;
	int localMaxima=numMaxima=0;
    
    //do a quick check to see if any one double homozygote has the highest fitness of all genotypes
    //if so, that corner is the solution
    int numBestGtypes=0;
    std::vector<std::vector<long double> > coordsOfBestGtype;
    for(int c=0;c<9;++c){std::vector<long double> coord(2,-one);coordsOfBestGtype.push_back(coord);}//error code
        

    long double wBarAtBestGtype=-one;

    std::vector<long double> wBarsToCompare(10,-one);//,x=-one,y=-one;-->make some extra slots
    std::vector<long double> pCoordsCompared(10,-one),qCoordsCompared(10,-one);
    wBarsToCompare[0]=wAABB;pCoordsCompared[0]=one;qCoordsCompared[0]=one;
    wBarsToCompare[1]=wAABb;pCoordsCompared[1]=one;qCoordsCompared[1]=half;
    wBarsToCompare[2]=wAAbb;pCoordsCompared[2]=one;qCoordsCompared[2]=zero;
    wBarsToCompare[3]=wAaBB;pCoordsCompared[3]=half;qCoordsCompared[3]=one;
    wBarsToCompare[4]=wAaBb;pCoordsCompared[4]=half;qCoordsCompared[4]=half;
    wBarsToCompare[5]=wAabb;pCoordsCompared[5]=half;qCoordsCompared[5]=zero;
    wBarsToCompare[6]=waaBB;pCoordsCompared[6]=zero;qCoordsCompared[6]=one;
    wBarsToCompare[7]=waaBb;pCoordsCompared[7]=zero;qCoordsCompared[7]=half;
    wBarsToCompare[8]=waabb;pCoordsCompared[8]=zero;qCoordsCompared[8]=zero;
    CoordinatesAtMax(wBarsToCompare,pCoordsCompared,qCoordsCompared,9,wBarAtBestGtype,coordsOfBestGtype,numBestGtypes);//compare 9
    if(numBestGtypes>1){//just in case
        CullDuplicateCoordinates2D(coordsOfBestGtype, numBestGtypes);}

    if(numBestGtypes==1){
            bool cornerSolution=false;
            if(coordsOfBestGtype[0][0]==one && coordsOfBestGtype[0][1]==one){
                    cornerSolution=true; wBarMax=wAABB;}
                else if(coordsOfBestGtype[0][0]==one && coordsOfBestGtype[0][1]==zero){
                    cornerSolution=true; wBarMax=wAAbb;}
                else if(coordsOfBestGtype[0][0]==zero && coordsOfBestGtype[0][1]==one){
                   cornerSolution=true; wBarMax=waaBB;}
                else if(coordsOfBestGtype[0][0]==zero && coordsOfBestGtype[0][1]==zero){
                    cornerSolution=true; wBarMax=waabb;}
                else{}
            if(cornerSolution==true){
                phat[0]=coordsOfBestGtype[0][0]; qhat[0]=coordsOfBestGtype[0][1];
                numMaxima=1;
                solutionFound=true;}
            }//numBestGtypes==1
        else if(numBestGtypes==9){//entirely flat fitness surface
            phat[0]=qhat[0]=-one;//return as a single solution instead of many
            pNeutral[0]=qNeutral[0]=true;
            numMaxima=1;
            solutionFound=true;
            }
        else{}
    
    if(solutionFound){
        return 1;}

    
    //define variables for values at the 4 corners {p,q} = {0,0}, {0,1}, {1,0}, & {1,1}
    long double pAt11=one,pAt01=zero,pAt10=one,pAt00=zero;
    long double qAt11=one,qAt01=one,qAt10=zero,qAt00=zero;
    long double pAthalf0=half,pAthalf1=half,pAt0half=zero,pAt1half=one;
    long double qAthalf0=zero,qAthalf1=one,qAt0half=half,qAt1half=half;
    long double wBarAt00,wBarAt01,wBarAt10,wBarAt11;
    long double wBarAthalf0,wBarAthalf1,wBarAt0half,wBarAt1half;//,wBarAtHalfHalf;
    long double dWbar_dpAt00,dWbar_dpAt01,dWbar_dpAt10,dWbar_dpAt11;
    long double dWbar_dqAt00,dWbar_dqAt01,dWbar_dqAt10,dWbar_dqAt11;
    long double dWbar_dpAt0half,dWbar_dpAt1half;
    long double dWbar_dqAthalf0,dWbar_dqAthalf1;
    bool pClimbsFromPt00,pClimbsFromPt01,pClimbsFromPt10,pClimbsFromPt11;//easier to interpret slopes
    bool qClimbsFromPt00,qClimbsFromPt01,qClimbsFromPt10,qClimbsFromPt11;
    bool pClimbsFromPt0half,pClimbsFromPt1half;
    bool qClimbsFromPthalf0,qClimbsFromPthalf1;

	wBarAt00=wBar(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt00,qAt00);
	dWbar_dpAt00=DwBar_dp(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt00,qAt00);
	dWbar_dqAt00=DwBar_dq(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt00,qAt00);
    dWbar_dpAt00=ROUND(dWbar_dpAt00,decimalDigitsToRound);
    dWbar_dqAt00=ROUND(dWbar_dqAt00,decimalDigitsToRound);
	pClimbsFromPt00 = (dWbar_dpAt00 > zero);//? use >= -> no
	qClimbsFromPt00 = (dWbar_dqAt00 > zero);

	wBarAt01=wBar(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt01,qAt01);
	dWbar_dpAt01=DwBar_dp(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt01,qAt01);
	dWbar_dqAt01=DwBar_dq(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt01,qAt01);
    dWbar_dpAt01=ROUND(dWbar_dpAt01,decimalDigitsToRound);
    dWbar_dqAt01=ROUND(dWbar_dqAt01,decimalDigitsToRound);
	pClimbsFromPt01 = (dWbar_dpAt01 > zero);
	qClimbsFromPt01 = (dWbar_dqAt01 < zero);// negative slope is up when q==1

	wBarAt10=wBar(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt10,qAt10);
	dWbar_dpAt10=DwBar_dp(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt10,qAt10);
	dWbar_dqAt10=DwBar_dq(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt10,qAt10);
    dWbar_dpAt10=ROUND(dWbar_dpAt10,decimalDigitsToRound);
    dWbar_dqAt10=ROUND(dWbar_dqAt10,decimalDigitsToRound);
	pClimbsFromPt10 = (dWbar_dpAt10 < zero);// negative slope is up when p==1
	qClimbsFromPt10 = (dWbar_dqAt10 > zero);

	wBarAt11=wBar(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt11,qAt11);
	dWbar_dpAt11=DwBar_dp(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt11,qAt11);
	dWbar_dqAt11=DwBar_dq(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt11,qAt11);
    dWbar_dpAt11=ROUND(dWbar_dpAt11,decimalDigitsToRound);
    dWbar_dqAt11=ROUND(dWbar_dqAt11,decimalDigitsToRound);
	pClimbsFromPt11 = (dWbar_dpAt11 < zero);// negative slope is up when p==1
	qClimbsFromPt11 = (dWbar_dqAt11 < zero);// negative slope is up when q==1

    wBarAthalf0=wBar(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAthalf0,qAthalf0);
    dWbar_dqAthalf0=DwBar_dq(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAthalf0,qAthalf0);
    dWbar_dqAthalf0=ROUND(dWbar_dqAthalf0,decimalDigitsToRound);
    qClimbsFromPthalf0 = (dWbar_dqAthalf0 > zero);//positive slope when p=0

    wBarAthalf1=wBar(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAthalf1,qAthalf1);
    dWbar_dqAthalf1=DwBar_dq(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAthalf1,qAthalf1);
    dWbar_dqAthalf1=ROUND(dWbar_dqAthalf1,decimalDigitsToRound);
    qClimbsFromPthalf1 = (dWbar_dqAthalf1 < zero);// negative slope is up when q==1

    wBarAt0half=wBar(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt0half,qAt0half);
    dWbar_dpAt0half=DwBar_dp(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt0half,qAt0half);
    dWbar_dpAt0half=ROUND(dWbar_dpAt0half,decimalDigitsToRound);
    pClimbsFromPt0half = (dWbar_dpAt0half > zero);// positive slope when q=0

    wBarAt1half=wBar(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt1half,qAt1half);
    dWbar_dpAt1half=DwBar_dp(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pAt1half,qAt1half);
    dWbar_dpAt1half=ROUND(dWbar_dpAt1half,decimalDigitsToRound);
    pClimbsFromPt1half = (dWbar_dpAt1half < zero);// negative slope is up when q==1

//    wBarAtHalfHalf=wBar(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,half,half);


	int numBestCorners=0;
	std::vector<std::vector<long double> > coordsAtBestCorners;
	for(int c=0;c<10;++c){std::vector<long double> coord(2,-one);coordsAtBestCorners.push_back(coord);}//error code

    long double wBarAtBestCorner=-one;

//	long double wBarsToCompare[10];//,x=-one,y=-one;-->make some extra slots
//	long double pCoordsCompared[10],qCoordsCompared[10];
	wBarsToCompare[0]=wBarAt00;pCoordsCompared[0]=pAt00;qCoordsCompared[0]=qAt00;
	wBarsToCompare[1]=wBarAt01;pCoordsCompared[1]=pAt01;qCoordsCompared[1]=qAt01;
	wBarsToCompare[2]=wBarAt10;pCoordsCompared[2]=pAt10;qCoordsCompared[2]=qAt10;
	wBarsToCompare[3]=wBarAt11;pCoordsCompared[3]=pAt11;qCoordsCompared[3]=qAt11;
	for(int i=4;i<10;++i){wBarsToCompare[i]=pCoordsCompared[i]=qCoordsCompared[i]=-one;}//initialize the rest with error codes
		
	CoordinatesAtMax(wBarsToCompare,pCoordsCompared,qCoordsCompared,4,wBarAtBestCorner,coordsAtBestCorners,numBestCorners);//compare just the first 4
	if(numBestCorners>1){//just in case
		CullDuplicateCoordinates2D(coordsAtBestCorners, numBestCorners);}

	

	//if all corners point downhill, then all single-het fitnesses are lower than their corners
		//and either (i) one or more corners hold wbarMax,
		//or (ii) there's a peak in the middle that may or may not be higher than the best corner
	if(!pClimbsFromPt00 && !qClimbsFromPt00 && !pClimbsFromPt01 && !qClimbsFromPt01 && !pClimbsFromPt10 && !qClimbsFromPt10 && !pClimbsFromPt11 && !qClimbsFromPt11){
		//handle case (i) here; if it's eliminated then handle case (ii) later
        if(!(pClimbsFromPt0half && pClimbsFromPt1half && qClimbsFromPthalf0 && qClimbsFromPthalf1)){
			//fitness doesn't climb from any edge toward the 2D center
			//so there can't be a central peak; one or more corners must have the highest fitness
			wBarMax=wBarAtBestCorner;
			numMaxima=numBestCorners;
			for(int c=0;c<numBestCorners;++c){
				phat[c]=coordsAtBestCorners[c][0]; qhat[c]=coordsAtBestCorners[c][1];}
			solutionFound=true;
			}//double-het fitness is lower than all 4 edges; completely convex landscape
		}//all corners pointed downhill
	if(!solutionFound){
		//if any corners point uphill along either edge, then the maximum fitness can be on that edge,
			//a different edge, or the middle
		//check for convex edges
		if((pClimbsFromPt00 && pClimbsFromPt10) || (qClimbsFromPt00 && qClimbsFromPt01)
				|| (pClimbsFromPt01 && pClimbsFromPt11) || (qClimbsFromPt10 && qClimbsFromPt11)){//there's at least one convex edge
				//this handles cases where the landscape slopes up to the middle of at least one edge,
				//as well as saddle cases, where two edges have higher fitness than the center such
				//that there are local maxima on more than one edge
			//there's a saddle if more than one edge is convex and the derivatives at their peaks all point down
				//if only one convex edge, then it will have the highest fitness
					//if wAaBb<fitness of that edge's heterozyote
			//check each edge to see if its het is better than the double-het
			long double wBarAtBestEdge=-one, wbarEdge=-one;
			std::vector<long double> pEdges(8,-one);
			std::vector<long double> qEdges(8,-one);
            std::vector<bool> pNeutralEdges(8,false);
            std::vector<bool> qNeutralEdges(8,false);
			std::vector<std::vector<long double> > coordsAtBestEdges;//4 edges, potentially maximized at both corners
			for(int edge=0;edge<8;++edge){std::vector<long double> coord(2,-one);coordsAtBestEdges.push_back(coord);}
			int edgeMaximaToConsider=0;
                    //this strategy doesn't work, because in some saddles, the edge maxima can be off center enough that
                        //the derivatives at those edges' midpoints might still point up.
                    //need to maximize each edge, apparently, and check that against wBar at {1/2,1/2}
            if((qClimbsFromPt00 && qClimbsFromPt01)){// && !pClimbsFromPt0half){//not 2-D convex from {0,1/2}
				//the p=0 edge might be maximal
				MaximizePopMeanFitnessQ(waaBB,waaBb,waabb,wbarEdge,qEdges,qNeutralEdges,localMaxima);
                long double dAtMax = DwBar_dp(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,zero,qEdges[0]);
                if(dAtMax<zero){
                    for(int m=0;m<localMaxima;++m){
                        wBarsToCompare[edgeMaximaToConsider]=wbarEdge;
                        pCoordsCompared[edgeMaximaToConsider]=zero;qCoordsCompared[edgeMaximaToConsider]=qEdges[m];
                        edgeMaximaToConsider++;}}
				}//p=0 edge
            if((qClimbsFromPt10 && qClimbsFromPt11)){// && !pClimbsFromPt1half){//not 2-D convex from {1,1/2}
				//the p=1 edge might be maximal
				MaximizePopMeanFitnessQ(wAABB,wAABb,wAAbb,wbarEdge,qEdges,qNeutralEdges,localMaxima);
                long double dAtMax = DwBar_dp(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,one,qEdges[0]);
                if(dAtMax>zero){//slopes down toward center
                    for(int m=0;m<localMaxima;++m){
                        wBarsToCompare[edgeMaximaToConsider]=wbarEdge;
                        pCoordsCompared[edgeMaximaToConsider]=one;qCoordsCompared[edgeMaximaToConsider]=qEdges[m];
                        edgeMaximaToConsider++;}}
				}//p=1 edge
            if((pClimbsFromPt00 && pClimbsFromPt10)){// && !qClimbsFromPthalf0){//not 2-D convex from {1/2,0}
				//the q=0 edge might be maximal
				MaximizePopMeanFitnessP(wAAbb,wAabb,waabb,wbarEdge,pEdges,pNeutralEdges,localMaxima);
                long double dAtMax = DwBar_dq(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pEdges[0],zero);
                if(dAtMax<zero){
                    for(int m=0;m<localMaxima;++m){
                        wBarsToCompare[edgeMaximaToConsider]=wbarEdge;
                        pCoordsCompared[edgeMaximaToConsider]=pEdges[m];qCoordsCompared[edgeMaximaToConsider]=zero;
                        edgeMaximaToConsider++;}}
				}//q=0 edge
            if((pClimbsFromPt01 && pClimbsFromPt11)){// && !qClimbsFromPthalf1){//not 2-D convex from {1/2,1}
				//the q=1 edge might be maximal
				MaximizePopMeanFitnessP(wAABB,wAaBB,waaBB,wbarEdge,pEdges,pNeutralEdges,localMaxima);
                long double dAtMax = DwBar_dq(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pEdges[0],one);
                if(dAtMax>zero){//slopes down towards center
                    for(int m=0;m<localMaxima;++m){
                        wBarsToCompare[edgeMaximaToConsider]=wbarEdge;
                        pCoordsCompared[edgeMaximaToConsider]=pEdges[m];qCoordsCompared[edgeMaximaToConsider]=one;
                        edgeMaximaToConsider++;}}
				}//q=1 edge
			if(edgeMaximaToConsider>0){
				//there's at least one edge that's not 2-D convex toward the center
				int numBestEdges=0;
				CoordinatesAtMax(wBarsToCompare,pCoordsCompared,qCoordsCompared,edgeMaximaToConsider,
					wBarAtBestEdge,coordsAtBestEdges,numBestEdges);
				//there can still be a higher corner

				if(wBarAtBestEdge>wBarAtBestCorner){
						wBarMax=wBarAtBestEdge;
						numMaxima=numBestEdges;
						for(int edge=0;edge<numBestEdges;++edge){
							phat[edge]=coordsAtBestEdges[edge][0];
							qhat[edge]=coordsAtBestEdges[edge][1];}
						solutionFound=true;
						}
					else if(wBarAtBestEdge<wBarAtBestCorner){
						wBarMax=wBarAtBestCorner;
						numMaxima=numBestCorners;
						for(int c=0;c<numBestCorners;++c){
							phat[c]=coordsAtBestCorners[c][0]; qhat[c]=coordsAtBestCorners[c][1];}
						solutionFound=true;
						}
					else{//neutral loci: corners and edges have equal maxima (which might be duplicates)
						coutLock.lock();
                        std::cout<<"Error in MaximizePopMeanFitnessPandQv2(): can't handle cases where wbar is maximal both in corner(s) and edge(s)"<<std::endl;
                        std::cout<<"                               for ref gtype="<<focalGtypeStr<<std::endl;
                        coutLock.unlock();
                        if(wAABB==wAABb && wAABB==wAAbb){//cis locus is neutral
                            qhat[0]=two;
                            MaximizePopMeanFitnessP(wAABB,wAaBB,waaBB,wbarEdge,pEdges,pNeutralEdges,localMaxima);
                            }
                        if(waaBB==waaBb && waaBB==waabb){//cis locus is neutral
                            qhat[0]=two;
                            MaximizePopMeanFitnessP(waaBB,waaBb,waabb,wbarEdge,pEdges,pNeutralEdges,localMaxima);
                            }
                        if(wAABB==wAaBB && wAABB==waaBB){//TF locus is neutral
                            phat[0]=two;
                            MaximizePopMeanFitnessQ(wAABB,wAABb,wAAbb,wbarEdge,qEdges,qNeutralEdges,localMaxima);
                            }
                        if(wAAbb==wAabb && wAABB==waabb){//TF locus is neutral
                            phat[0]=two;
                            MaximizePopMeanFitnessQ(waaBB,waaBb,waabb,wbarEdge,qEdges,qNeutralEdges,localMaxima);
                            }
                        solutionFound=true;
						}//neutral loci
				}//edgeMaximaToConsider>0
			//there must be a peak in the middle of the distribution someplace; check that next
			}//there's a convex edge
		}// !solutionFound

    //maybe there are neutral edges with maximal fitness
    if(!solutionFound){
        if(numBestCorners==2){
                if(((coordsAtBestCorners[0][0]==0 && coordsAtBestCorners[0][1]==0)//pt {0,0}
                    && (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==1))//pt {0,1}
                    ||
                   ((coordsAtBestCorners[0][0]==0 && coordsAtBestCorners[0][1]==1)//pt {0,1}
                    && (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==0))){//pt {0,0}
                        if(!pClimbsFromPt0half){
                            if(waaBB==waaBb){//neutral q at p=0
                                    wBarMax=waaBB;
                                    qNeutral[0]=true;
                                    phat[0]=zero;
                                    qhat[0]=-one;//flag
                                    numMaxima=1;
                                    solutionFound=true;}
                                else{//both corners are solutions
                                    wBarMax=waabb;
                                    phat[0]=zero; qhat[0]=zero;
                                    phat[1]=zero; qhat[1]=one;
                                    numMaxima=2;
                                    solutionFound=true;}
                            }//waaBb>=wAaBb
                        }// {0,0} & {0,1}
                    else if(((coordsAtBestCorners[0][0]==0 && coordsAtBestCorners[0][1]==1)//pt {0,1}
                        && (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==1))//pt {1,1}
                        ||
                        ((coordsAtBestCorners[0][0]==1 && coordsAtBestCorners[0][1]==1)//pt {1,1}
                        && (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==1))){//pt {0,1}
                        if(!qClimbsFromPthalf1){
                            if(wAABB==wAaBB){//neutral p at q=1
                                    wBarMax=wAABB;
                                    pNeutral[0]=true;
                                    qhat[0]=one;
                                    phat[0]=-one;//flag
                                    numMaxima=1;
                                    solutionFound=true;}
                                else{//both corners are solutions
                                    wBarMax=wAABB;
                                    phat[0]=zero; qhat[0]=one;
                                    phat[1]=one; qhat[1]=one;
                                    numMaxima=2;
                                    solutionFound=true;}
                            }//wAaBB>=wAaBb
                        }//{0,1} & {1,1}
                    else if(((coordsAtBestCorners[0][0]==1 && coordsAtBestCorners[0][1]==1)//pt {1,1}
                        && (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==0))//pt {1,0}
                        ||
                            ((coordsAtBestCorners[0][0]==1 && coordsAtBestCorners[0][1]==0)//pt {1,0}
                        && (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==1))){//pt {1,1}
                        if(!pClimbsFromPt1half){
                            if(wAABB==wAABb){//neutral q at p=1
                                    wBarMax=wAABB;
                                    qNeutral[0]=true;
                                    qhat[0]=-one;//flag
                                    phat[0]=one;
                                    numMaxima=1;
                                    solutionFound=true;}
                                else{//both corners are solutions
                                    wBarMax=wAABB;
                                    phat[0]=one; qhat[0]=zero;
                                    phat[1]=one; qhat[1]=one;
                                    numMaxima=2;
                                    solutionFound=true;}
                            }//wAABb>=wAaBb
                        }//{1,0} & {1,1}
                    else if(((coordsAtBestCorners[0][0]==1 && coordsAtBestCorners[0][1]==0)//pt {1,0}
                        && (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==0))//pt {0,0}
                        ||
                            ((coordsAtBestCorners[0][0]==0 && coordsAtBestCorners[0][1]==0)//pt {0,0}
                        && (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==0))){//pt {1,0}
                        if(!qClimbsFromPthalf0){
                            if(wAAbb==wAabb){//neutral p at q=0
                                    wBarMax=wAAbb;
                                    pNeutral[0]=true;
                                    qhat[0]=zero;
                                    phat[0]=-one;//flag
                                    numMaxima=1;
                                    solutionFound=true;}
                                else{//both corners are solutions
                                    wBarMax=wAAbb;
                                    phat[0]=zero; qhat[0]=zero;
                                    phat[1]=one; qhat[1]=zero;
                                    numMaxima=2;
                                    solutionFound=true;}
                            }//wAabb>=wAaBb
                        }//{0,0} & {1,0}
                    else{
                        if(((coordsAtBestCorners[0][0]==0 && coordsAtBestCorners[0][1]==0)//{0,0}
                           && (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==1))//{1,1}
                           ||
                           ((coordsAtBestCorners[0][0]==1 && coordsAtBestCorners[0][1]==1)//{1,1}
                           && (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==0))){//{0,0}
                                if(wAABB>=wAaBb){
                                    wBarMax=wAABB;
                                    phat[0]=zero; qhat[0]=zero;
                                    phat[1]=one; qhat[1]=one;
                                    numMaxima=2;
                                    solutionFound=true;}
                                }
                            else if(((coordsAtBestCorners[0][0]==0 && coordsAtBestCorners[0][1]==1)//{0,1}
                            && (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==0))//{1,0}
                            ||
                            ((coordsAtBestCorners[0][0]==1 && coordsAtBestCorners[0][1]==0)//{1,0}
                            && (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==1))){//{0,1}
                                if(wAAbb>=wAaBb){
                                    wBarMax=wAAbb;
                                    phat[0]=zero; qhat[0]=one;
                                    phat[1]=one; qhat[1]=zero;
                                    numMaxima=2;
                                    solutionFound=true;}
                                }
                            else{
                                coutLock.lock(); std::cout<<"shouldn't reach this point"<<std::endl; coutLock.unlock();}
                        }//equal & not adjacent
                        
                    }//numBestCorners==2
            else if(numBestCorners==3){
                //first identify the corners, then check the edges
                //{0,0},{0,1},{1,1}
                if(   ((coordsAtBestCorners[0][0]==0 && coordsAtBestCorners[0][1]==0)//{0,0}
                           &&                           // and either
                            ( ( (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==1)
                                    &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==1)//( {0,1} & {1,1}
                                )
                             ||
                               ( (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==1)//or
                                    &&(coordsAtBestCorners[2][0]==0 && coordsAtBestCorners[2][1]==1))//{1,1} & {0,1} )
                              ))
                        || //or a different sequence
                           ((coordsAtBestCorners[0][0]==0 && coordsAtBestCorners[0][1]==1)//{0,1}
                           &&                           // and either
                            ( ( (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==1)
                                    &&(coordsAtBestCorners[2][0]==0 && coordsAtBestCorners[2][1]==0)//( {1,1} & {0,0}
                                )
                             ||
                               ( (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==0)//or
                                    &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==1))//{0,0} & {1,1} )
                              ))
                       || //or a different sequence
                          ((coordsAtBestCorners[0][0]==1 && coordsAtBestCorners[0][1]==1)//{1,1}
                          &&                           // and either
                           ( ( (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==1)
                                   &&(coordsAtBestCorners[2][0]==0 && coordsAtBestCorners[2][1]==0)//( {0,1} & {0,0}
                               )
                            ||
                              ( (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==0)//or
                                   &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==1))//{0,0} & {0,1} )
                             ))
                        ){//corners are aabb, aaBB and AABB
                            //potential cases:
                                //just the corners are high
                                //one or more edges is neutral
                                //{0,0},{0,1},{1,1}
                        if(waabb>=wAaBb){
                            if(waabb == waaBb && wAABB == wAaBB){//both edges neutral when the adjacent edge is fixed
                                    wBarMax=waaBB;
                                    qNeutral[0]=true;//q is neutral when p=0
                                    phat[0]=zero;
                                    qhat[0]=-one;
                                    pNeutral[1]=true;//p is neutral when q=1
                                    qhat[1]=one;
                                    phat[1]=-one;
                                    localMaxima=2;
                                    }
                                else if(waabb == waaBb){//q is neutral when p=0 and the {1,1} corner is a solution
                                    wBarMax=wAABB;
                                    phat[0]=one;
                                    qhat[0]=one;
                                    qNeutral[1]=true;
                                    phat[1]=zero;
                                    qhat[1]=-one;
                                    localMaxima=2;}
                                else if(wAABB == wAaBB){//p is neutral when q=1 and the {0,0} corner is a solution
                                    wBarMax=waabb;
                                    phat[0]=zero;
                                    qhat[0]=zero;
                                    pNeutral[1]=true;
                                    phat[1]=-one;
                                    qhat[1]=one;
                                    localMaxima=2;}
                                else{//all 3 corners are unique solutions
                                    wBarMax=waabb;
                                    phat[0]=zero;
                                    qhat[0]=zero;
                                    phat[1]=zero;
                                    qhat[1]=one;
                                    phat[2]=one;
                                    qhat[2]=one;
                                    localMaxima=3;}
                            solutionFound=true;
                            }
                        }
                //{0,1},{1,1},(1,0}
                else if( ((coordsAtBestCorners[0][0]==0 && coordsAtBestCorners[0][1]==1)//{0,1}
                        &&                           // and either
                         ( ( (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==1)
                                 &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==0)//( {1,1} & {1,0}
                             )
                          ||
                            ( (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==0)//or
                                 &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==1))//{1,0} & {1,1} )
                           ))
                     || //or a different sequence
                        ((coordsAtBestCorners[0][0]==1 && coordsAtBestCorners[0][1]==1)//{1,1}
                        &&                           // and either
                         ( ( (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==1)
                                 &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==0)//( {0,1} & {1,0}
                             )
                          ||
                            ( (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==0)//or
                                 &&(coordsAtBestCorners[2][0]==0 && coordsAtBestCorners[2][1]==1))//{1,0} & {0,1} )
                           ))
                    || //or a different sequence
                       ((coordsAtBestCorners[0][0]==1 && coordsAtBestCorners[0][1]==0)//{1,0}
                       &&                           // and either
                        ( ( (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==1)
                                &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==1)//( {0,1} & {1,1}
                            )
                         ||
                           ( (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==1)//or
                                &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==1))//{1,1} & {0,1} )
                          ))
                     ){//corners are aaBB, AABB and AAbb
                         //potential cases:
                             //just the corners are high
                             //one or more edges is neutral
                        //{0,1},{1,1},(1,0}
                        if(waaBB>=wAaBb){
                            if(waaBB == wAaBB && wAAbb == wAABb){//both edges neutral when the adjacent edge is fixed
                                    wBarMax=wAAbb;
                                    pNeutral[0]=true;//p is neutral when q=0
                                    phat[0]=-one;
                                    qhat[0]=zero;
                                    qNeutral[1]=true;//q is neutral when p=1
                                    phat[1]=one;
                                    qhat[1]=-one;
                                    localMaxima=2;
                                    }
                                else if(waaBB == wAaBB){//p is neutral when q=0, and the {1,0} corner is a solution
                                    wBarMax=wAAbb;
                                    phat[0]=one;
                                    qhat[0]=zero;
                                    pNeutral[1]=true;
                                    phat[1]=-one;
                                    qhat[1]=zero;
                                    localMaxima=2;
                                    }
                                else if(wAAbb == wAABb){//q is neutral when p=1, and the {0,1} corner is a solution
                                    wBarMax=waaBB;
                                    phat[0]=zero;
                                    qhat[0]=one;
                                    qNeutral[1]=true;
                                    phat[1]=one;
                                    qhat[1]=-one;
                                    localMaxima=2;
                                    }
                                else{//all 3 corners are unique solutions
                                    wBarMax=waaBB;
                                    phat[0]=zero;
                                    qhat[0]=one;
                                    phat[1]=one;
                                    qhat[1]=zero;
                                    phat[2]=one;
                                    qhat[2]=one;
                                    localMaxima=3;}
                            solutionFound=true;
                            }
                        }
                //{1,1},{1,0},{0,0}
                else if( ((coordsAtBestCorners[0][0]==1 && coordsAtBestCorners[0][1]==1)//{1,1}
                           &&                           // and either
                            ( ( (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==0)
                                    &&(coordsAtBestCorners[2][0]==0 && coordsAtBestCorners[2][1]==0)//( {1,0} & {0,0}
                                )
                             ||
                               ( (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==0)//or
                                    &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==0))//{0,0} & {1,0} )
                              ))
                        || //or a different sequence
                           ((coordsAtBestCorners[0][0]==1 && coordsAtBestCorners[0][1]==0)//{1,0}
                           &&                           // and either
                            ( ( (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==1)
                                    &&(coordsAtBestCorners[2][0]==0 && coordsAtBestCorners[2][1]==0)//( {1,1} & {0,0}
                                )
                             ||
                               ( (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==0)//or
                                    &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==1))//{0,0} & {1,1} )
                              ))
                       || //or a different sequence
                          ((coordsAtBestCorners[0][0]==0 && coordsAtBestCorners[0][1]==0)//{0,0}
                          &&                           // and either
                           ( ( (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==1)
                                   &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==0)//( {1,1} & {1,0}
                               )
                            ||
                              ( (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==0)//or
                                   &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==1))//{1,0} & {1,1} )
                             ))
                        ){//corners are AABB, AAbb and aabb
                            //{0,0},{1,0},{1,1}
                            //potential cases:
                                //just the corners are high
                                //one or more edges is neutral
                        if(waabb>=wAaBb){
                            if(waabb == wAabb && wAABB == wAABb){//both edges neutral when the adjacent edge is fixed
                                    wBarMax=wAAbb;
                                    qNeutral[0]=true;//q is neutral when p=1
                                    phat[0]=one;
                                    qhat[0]=-one;
                                    pNeutral[1]=true;//p is neutral when q=0
                                    phat[1]=-one;
                                    qhat[1]=zero;
                                    localMaxima=2;
                                    }
                                else if(wAABB == wAABb){//q is neutral when p=1
                                    //and the {0,0} corner is a solution
                                    wBarMax=waabb;
                                    phat[0]=zero;
                                    qhat[0]=zero;
                                    qNeutral[1]=true;
                                    phat[1]=one;
                                    qhat[1]=-one;
                                    localMaxima=2;
                                    }
                                else if(waabb == wAabb){//p is neutral when q=0
                                    //and the {1,1} corner is a solution
                                    wBarMax=wAABB;
                                    phat[0]=one;
                                    qhat[0]=one;
                                    pNeutral[1]=true;
                                    phat[1]=-one;
                                    qhat[1]=zero;
                                    localMaxima=2;
                                    }
                                else{//all 3 corners are unique solutions
                                    wBarMax=waabb;
                                    phat[0]=zero;
                                    qhat[0]=zero;
                                    phat[1]=one;
                                    qhat[1]=zero;
                                    phat[2]=one;
                                    qhat[2]=one;
                                    localMaxima=3;
                                    }
                            solutionFound=true;
                            }
                        }
                //{1,0},{0,0},{0,1}
                else if( ((coordsAtBestCorners[0][0]==1 && coordsAtBestCorners[0][1]==0)//{1,0}
                          &&                           // and either
                           ( ( (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==0)
                                   &&(coordsAtBestCorners[2][0]==0 && coordsAtBestCorners[2][1]==1)//( {0,0} & {0,1}
                               )
                            ||
                              ( (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==1)//or
                                   &&(coordsAtBestCorners[2][0]==0 && coordsAtBestCorners[2][1]==0))//{0,1} & {0,0} )
                             ))
                       || //or a different sequence
                          ((coordsAtBestCorners[0][0]==0 && coordsAtBestCorners[0][1]==0)//{0,0}
                          &&                           // and either
                           ( ( (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==0)
                                   &&(coordsAtBestCorners[2][0]==0 && coordsAtBestCorners[2][1]==1)//( {1,0} & {0,1}
                               )
                            ||
                              ( (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==1)//or
                                   &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==0))//{0,1} & {1,0} )
                             ))
                      || //or a different sequence
                         ((coordsAtBestCorners[0][0]==0 && coordsAtBestCorners[0][1]==1)//{0,1}
                         &&                           // and either
                          ( ( (coordsAtBestCorners[1][0]==1 && coordsAtBestCorners[1][1]==0)
                                  &&(coordsAtBestCorners[2][0]==0 && coordsAtBestCorners[2][1]==0)//( {1,0} & {0,0}
                              )
                           ||
                             ( (coordsAtBestCorners[1][0]==0 && coordsAtBestCorners[1][1]==0)//or
                                  &&(coordsAtBestCorners[2][0]==1 && coordsAtBestCorners[2][1]==0))//{0,0} & {1,0} )
                            ))
                       ){//corners are AAbb, aabb and aaBB
                           //{1,0},{0,0},{0,1}
                            //potential cases:
                               //just the corners are high
                               //one or more edges is neutral
                        if(wAAbb>=wAaBb){
                            if(wAAbb == wAabb && waaBB == waaBb){//both edges neutral when the adjacent edge is fixed
                                wBarMax=wAAbb;
                                qNeutral[0]=true;//q is neutral when p=1
                                phat[0]=one;
                                qhat[0]=-one;
                                pNeutral[1]=true;//p is neutral when q=0
                                phat[1]=-one;
                                qhat[1]=zero;
                                localMaxima=2;}
                           else if(wAAbb == wAabb){//q is neutral when p=1
                                //and the {0,1} corner is a solution
                                wBarMax=waaBB;
                                phat[0]=zero;
                                qhat[0]=one;
                                qNeutral[1]=true;
                                phat[1]=one;
                                qhat[1]=-one;
                                localMaxima=2;}
                           else if(waaBB == waaBb){//p is neutral when q=0
                                //and the {1,0} corner is a solution
                                wBarMax=wAAbb;
                                phat[0]=one;
                                qhat[0]=zero;
                                pNeutral[1]=true;
                                phat[1]=-one;
                                qhat[1]=zero;
                                localMaxima=2;}
                           else{//all 3 corners are unique solutions
                               wBarMax=waabb;
                               phat[0]=zero;
                               qhat[0]=zero;
                               phat[1]=zero;
                               qhat[1]=one;
                               phat[2]=one;
                               qhat[2]=zero;
                               localMaxima=3;}
                       solutionFound=true;
                       }
                    }
                else{
                    coutLock.lock();
                    std::cout<<"3-site case indicated but no 3-site combination found; shouldn't get here"<<std::endl;
                    coutLock.unlock();
                }
                //potential cases:
                    //just the corners are high
                    //one or more edges is neutral
                //more than one locus can be neutral when the other is not, but both can't be neutral sinultaneously
                }
            else if(numBestCorners==4){
                coutLock.lock(); std::cout<<"4 best corners"<<std::endl; coutLock.unlock();
                //potential cases:
                
                    if(waabb>=wAaBb){
                        if(waabb>wAABb && waabb>wAaBB && waabb>wAabb && waabb>waaBb){//just the corners are high
                            phat[0]=zero; qhat[0]=zero;
                            phat[1]=one; qhat[1]=zero;
                            phat[2]=zero; qhat[2]=one;
                            phat[3]=one; qhat[3]=one;
                            localMaxima=4;
                            solutionFound=true;
                            }
                        else if(waabb==wAABb && waabb>wAaBB && waabb>wAabb && waabb>waaBb){//1 neutral edge
                            //q is neutral edge at p=1; the two p=0 corners are solutions
                            wBarMax=waabb;
                            phat[0]=zero; qhat[0]=zero;
                            phat[1]=zero; qhat[1]=one;
                            phat[2]=one; qhat[2]=-one; qNeutral[2]=true;
                            localMaxima=3;
                            solutionFound=true;
                            }
                        else if(waabb>wAABb && waabb==wAaBB && waabb>wAabb && waabb>waaBb){//1 neutral edge
                              //p is neutral edge at q=1; the two q=0 corners are solutions
                              wBarMax=waabb;
                              phat[0]=zero; qhat[0]=zero;
                              phat[1]=one; qhat[1]=zero;
                              phat[2]=-one; qhat[2]=one; pNeutral[2]=true;
                              localMaxima=3;
                            solutionFound=true;
                              }
                        else if(waabb>wAABb && waabb>wAaBB && waabb==wAabb && waabb>waaBb){//1 neutral edge
                            //p is neutral edge at q=0; the two q=1 corners are solutions
                            wBarMax=waaBB;
                            phat[0]=zero; qhat[0]=one;
                            phat[1]=one; qhat[1]=one;
                            phat[2]=-one; qhat[2]=zero; pNeutral[2]=true;
                            localMaxima=3;
                            solutionFound=true;
                            }
                        else if(waabb>wAABb && waabb>wAaBB && waabb>wAabb && waabb==waaBb){//1 neutral edge
                            //q is neutral edge at p=0; the two p=1 corners are solutions
                            wBarMax=wAAbb;
                            phat[0]=one; qhat[0]=zero;
                            phat[1]=one; qhat[1]=one;
                            phat[2]=zero; qhat[2]=-one; qNeutral[2]=true;
                            localMaxima=3;
                            solutionFound=true;
                            }
                        else if(waabb==wAABb && waabb==wAaBB && waabb>wAabb && waabb>waaBb){//2 neutral edges
                            //q is neutral edge at p=1 ({1,0} & {1,1}); p is neutral edge at q=1 ({0,1} & {1,1}); the {0,0} corner is a solution
                            wBarMax=waabb;
                            phat[0]=zero; qhat[0]=zero;
                            phat[1]=one; qhat[1]=-one; qNeutral[1]=true;
                            phat[2]=-one; qhat[2]=one; pNeutral[2]=true;
                            localMaxima=3;
                            solutionFound=true;
                            }
                        else if(waabb==wAABb && waabb>wAaBB && waabb==wAabb && waabb>waaBb){//2 neutral edges
                            //q is neutral edge at p=1 ({1,0} & {1,1}); p is neutral edge at q=0 ({0,0} & {1,0}); the {0,1} corner is a solution
                            wBarMax=waaBB;
                            phat[0]=zero; qhat[0]=one;
                            phat[1]=one; qhat[1]=-one; qNeutral[1]=true;
                            phat[2]=-one; qhat[2]=zero; pNeutral[2]=true;
                            localMaxima=3;
                            solutionFound=true;
                            }
                        else if(waabb==wAABb && waabb>wAaBB && waabb>wAabb && waabb==waaBb){//2 neutral edges
                            //q is neutral edge at p=1 ({1,0} & {1,1}); q is neutral edge at p=0 ({0,0} & {0,1}); there are no corner solutions
                            wBarMax=waaBB;
                            phat[0]=zero; qhat[0]=-one; qNeutral[0]=true;
                            phat[1]=one; qhat[1]=-one; qNeutral[1]=true;
                            localMaxima=2;
                            solutionFound=true;
                            }
                       else if(waabb>wAABb && waabb==wAaBB && waabb==wAabb && waabb>waaBb){//2 neutral edges
                           //p is neutral edge at q=1 ({0,1} & {1,1}); p is neutral edge at q=0 ({0,0} & {1,0}); there are no corner solutions
                           wBarMax=wAAbb;
                           phat[0]=-one; qhat[0]=zero; pNeutral[0]=true;
                           phat[1]=-one; qhat[1]=one; pNeutral[1]=true;
                           localMaxima=2;
                           solutionFound=true;
                           }
                        else if(waabb>wAABb && waabb==wAaBB && waabb>wAabb && waabb==waaBb){//2 neutral edges
                            //p is neutral edge at q=1 ({0,1} & {1,1}); q is neutral edge at p=0 ({0,0} & {0,1}); the {1,0} corner is a solution
                            wBarMax=wAAbb;
                            phat[0]=one; qhat[0]=zero;
                            phat[1]=-one; qhat[1]=one; pNeutral[1]=true;
                            phat[2]=zero; qhat[2]=-one; qNeutral[2]=true;
                            localMaxima=3;
                            solutionFound=true;
                            }
                        else if(waabb>wAABb && waabb>wAaBB && waabb==wAabb && waabb==waaBb){//2 neutral edges
                            //p is neutral edge at q=0 ({0,0} & {1,0}) q is neutral edge at p=0 ({0,0} & {0,1}); the {1,1} corner is a solution
                            wBarMax=wAABB;
                            phat[0]=one; qhat[0]=one;
                            phat[1]=-one; qhat[1]=zero; pNeutral[1]=true;
                            phat[2]=zero; qhat[2]=-one; qNeutral[2]=true;
                            localMaxima=3;
                            solutionFound=true;
                            }
                        else if(waabb==wAABb && waabb==wAaBB && waabb==wAabb && waabb>waaBb){//3 neutral edges
                            //q is neutral edge at p=1 ({1,0} & {1,1})
                            //p is neutral edge at q=1 ({0,1} & {1,1})
                            //p is neutral edge at q=0 ({0,0} & {0,1})
                            wBarMax=wAABB;
                            phat[0]=one; qhat[0]=-one; qNeutral[0]=true;
                            phat[1]=-one; qhat[1]=zero; pNeutral[1]=true;
                            phat[2]=-one; qhat[2]=one; pNeutral[2]=true;
                            localMaxima=3;
                            solutionFound=true;
                            }
                        else if(waabb==wAABb && waabb==wAaBB && waabb>wAabb && waabb==waaBb){//3 neutral edges
                            //q is neutral edge at p=1 ({1,0} & {1,1})
                            //p is neutral edge at q=1 ({0,1} & {1,1})
                            //q is neutral edge at p=0 ({0,0} & {0,1})
                            wBarMax=waabb;
                            phat[0]=zero; qhat[0]=-one; qNeutral[0]=true;
                            phat[1]=one; qhat[1]=-one; qNeutral[1]=true;
                            phat[2]=-one; qhat[2]=one; pNeutral[2]=true;
                            localMaxima=3;
                            solutionFound=true;
                            }
                        else if(waabb==wAABb && waabb>wAaBB && waabb==wAabb && waabb==waaBb){//3 neutral edges
                            //q is neutral edge at p=1 ({1,0} & {1,1})
                            //p is neutral edge at q=0 ({0,1} & {1,1})
                            //q is neutral edge at p=0 ({0,0} & {0,1})
                            wBarMax=wAABB;
                            phat[0]=zero; qhat[0]=-one; pNeutral[0]=true;
                            phat[1]=one; qhat[1]=-one; qNeutral[1]=true;
                            phat[2]=-one; qhat[2]=zero; pNeutral[2]=true;
                            localMaxima=3;
                            solutionFound=true;
                            }
                        else if(waabb>wAABb && waabb==wAaBB && waabb==wAabb && waabb==waaBb){//3 neutral edges
                            //p is neutral edge at q=1 ({0,1} & {1,1})
                            //p is neutral edge at q=0 ({0,1} & {1,1})
                            //q is neutral edge at p=0 ({0,0} & {0,1})
                            wBarMax=wAABB;
                            phat[0]=zero; qhat[0]=-one; qNeutral[0]=true;
                            phat[1]=-one; qhat[1]=zero; pNeutral[1]=true;
                            phat[2]=-one; qhat[2]=one; pNeutral[2]=true;
                            localMaxima=3;
                            solutionFound=true;
                            }
                        else if(waabb==wAABb && waabb==wAaBB && waabb==wAabb && waabb==waaBb){//4 neutral edges
                            //q is neutral edge at p=1 ({1,0} & {1,1})
                            //p is neutral edge at q=1 ({0,1} & {1,1})
                            //p is neutral edge at q=0 ({0,1} & {1,1})
                            //q is neutral edge at p=0 ({0,0} & {0,1})
                            wBarMax=wAABB;
                            phat[0]=zero; qhat[0]=-one; qNeutral[0]=true;
                            phat[1]=one; qhat[1]=-one; qNeutral[1]=true;
                            phat[2]=-one; qhat[2]=zero; pNeutral[2]=true;
                            phat[3]=-one; qhat[3]=one; pNeutral[3]=true;
                            localMaxima=4;
                            solutionFound=true;
                            }
                        else{}
                    }
                }//4 best corners
            else{}
        }// !solutionFound
    
	if(!solutionFound){
		//the maximum might be in a local peak in the middle if derivatives from edges point up
			//but even so, a convex center can still be lower than the best corner
		//so, maximize the middle of the landscape, then compare it to the best corner

		std::vector<long double> coordsA(2,half), coordsB(2,-one), coordsC(2,-one);//{p,q} for these points
		int cornerMinima=0;
		std::vector<std::vector<long double> > coordsAtWorstCorner;//{p,q} for these points
		for(int c=0;c<4;++c){std::vector<long double> coord(2,-one);coordsAtWorstCorner.push_back(coord);}
		long double wBarAtWorstCorner=-one;
		CoordinatesAtMin(wBarsToCompare,pCoordsCompared,qCoordsCompared,4,wBarAtWorstCorner,coordsAtWorstCorner,cornerMinima);
		if(cornerMinima>1){//just in case
			CullDuplicateCoordinates2D(coordsAtWorstCorner,cornerMinima);}

		int nFuntionCalls=0;
		long double convergenceBracket=(long double)10.0e-6;
		long double wbarABC[3];
		long double pABC[3];
		long double qABC[3];
		//initialize at arbitrary points near the center
		pABC[0]=(long double)0.42;qABC[0]=(long double)0.45;//point A
		pABC[1]=(long double)0.41;qABC[1]=(long double)0.52;//point B
		pABC[2]=(long double)0.55; qABC[2]=(long double)0.55;//point C
		wbarABC[0] = wBar(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pABC[0],qABC[0]);
		wbarABC[1] = wBar(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pABC[1],qABC[1]);
		wbarABC[2] = wBar(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,pABC[2],qABC[2]);
		SortByWbar(wbarABC,pABC,qABC);
		int converged=amoeba(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,
									 wbarABC,pABC,qABC,convergenceBracket,nFuntionCalls);
		if(!converged){
			coutLock.lock(); std::cout<<"MaximizePopMeanFitnessPandQv2(...) failed to converge"<<std::endl; coutLock.unlock();}
		long double wBarAtBestMiddle,coordsAtHighestMiddle[2];
		wBarAtBestMiddle=wbarABC[0];
		coordsAtHighestMiddle[0]=pABC[0];
		coordsAtHighestMiddle[1]=qABC[0];
		wBarAtBestMiddle=MIN(MAX(ROUND(wBarAtBestMiddle,decimalDigitsToRound),zero),one);//round & bound
		coordsAtHighestMiddle[0]=MIN(MAX(ROUND(coordsAtHighestMiddle[0],decimalDigitsToRound),zero),one);
		coordsAtHighestMiddle[1]=MIN(MAX(ROUND(coordsAtHighestMiddle[1],decimalDigitsToRound),zero),one);
		//it may be that the center is lower that the best corner
        //or the maximization reached a high, neutral edge
        
        //see if the maximum reached an edge
        if(coordsAtHighestMiddle[0]==zero){//see if q is neutral on the p=0 edge
                if(waaBB==waaBb && waaBB==waabb && !pClimbsFromPt0half){//flat edge & the center isn't higher
                    wBarMax=wAAbb;
                    numMaxima=1;
                    phat[0]=zero;
                    qhat[0]=-one;
                    qNeutral[0]=true;}
                }//p==0
            else if(coordsAtHighestMiddle[0]==one){//see if q is neutral on the p=1 edge
                if(wAABB==wAABb && wAABB==wAAbb && !pClimbsFromPt1half){//flat edge & the center isn't higher
                    wBarMax=wAABB;
                    numMaxima=1;
                    phat[0]=one;
                    qhat[0]=-one;
                    qNeutral[0]=true;}
                }//p==1
            else if(coordsAtHighestMiddle[1]==zero){//see if p is neutral on the q=0 edge
                if(wAAbb==wAabb && wAAbb==waabb && !qClimbsFromPthalf0){//flat edge & the center isn't higher
                    wBarMax=wAAbb;
                    numMaxima=1;
                    phat[0]=-one;
                    qhat[0]=zero;
                    pNeutral[0]=true;}
                }//q==0
            else if(coordsAtHighestMiddle[1]==one){//see if p is neutral on the q=1 edge
                if(wAABB==wAaBB && wAABB==waaBB && !qClimbsFromPthalf1){//flat edge & the center isn't higher
                    wBarMax=wAABB;
                    numMaxima=1;
                    phat[0]=-one;
                    qhat[0]=one;
                    pNeutral[0]=true;}
                }//q==1
            else{
                //check if center maximum < corner maximum
                for(int i=0;i<numBestCorners;++i){
                    wBarsToCompare[i]=MIN(MAX(ROUND(wBarAtBestCorner,decimalDigitsToRound),zero),one);
                    pCoordsCompared[i]=coordsAtBestCorners[i][0];
                    qCoordsCompared[i]=coordsAtBestCorners[i][1];
                    }
                wBarsToCompare[numBestCorners]=wBarAtBestMiddle;
                pCoordsCompared[numBestCorners]=coordsAtHighestMiddle[0];
                qCoordsCompared[numBestCorners]=coordsAtHighestMiddle[1];
                //set wBarMAx, p and q
                int localMaxima=0, coordsToCompare=numBestCorners+1;
                CoordinatesAtMax(wBarsToCompare,pCoordsCompared,qCoordsCompared,coordsToCompare,wBarMax,
                    coordsAtBestCorners,localMaxima);
                if(localMaxima>1){//just in case
                    int oldLocalMaxima=localMaxima;
                    CullDuplicateCoordinates2D(coordsAtBestCorners,localMaxima);
                    for(int m=localMaxima;m<oldLocalMaxima;++m){
                        phat[m]=qhat[m]=-one;}//reset these since -one is used as an error flag
                    }
                for(int m=0;m<localMaxima;++m){
                    phat[m]=coordsAtBestCorners[m][0];
                    qhat[m]=coordsAtBestCorners[m][1];}
                numMaxima=localMaxima;
                }//center maximum ?< corner maximum
		}//!solutionFound

	return 1;
	}//MaximizePopMeanFitnessPandQv2()



void MaximizeUsingBitstringsOneReferenceGtype(SimplestRegPathIndividual* focalIndivP, simulationSettings* simSetP,
        genotypeSettings* gtypeSetP, FitnessMaximumSolutionSet* fmssP, long double maxPopMeanFitness,
        int* solutionEqualsMaxOrBetterP){

//int bitstringLen,
//            long double NtfsatPerAllele, long double deltaG1dosage, long double deltaG1, long double minExpression,
//            long double maxExpression, long double Popt, long double omega)
    //this takes a single 2-locus, 2-allele reference genotype, partitions it into genotypes comprising all possible
    //combinations of homozygous and heterozygous allele combinations at the 2 loci,
    //and calculates the equilibrium allele frequencies.  If the populaiton mean fitness is a new maximum,
    //it replaces updates the maximal fitness and stores the reference; if it equals the old maximun then it updates the
    //counters of the number of maximal-fitness reference genotypes; otherwise it records nothing.
    
    //it's constructed this way so that it can be handled in a single thrad for parallel computing

    //SimplestRegPathIndividual focalIndiv(false);//indivAaBb
    //focalIndiv.SetGenotype(2,1,cisVal1);//set before entering this function

    //passing pointers to this function gets around a compiler problem; dereference them here for readability
    SimplestRegPathIndividual& focalIndiv = *focalIndivP;
    simulationSettings& simSet = *simSetP;
    genotypeSettings& gtypeSet = *gtypeSetP;
    FitnessMaximumSolutionSet& fmss = *fmssP;
    int& solutionEqualsMaxOrBetter = *solutionEqualsMaxOrBetterP;
    bool splitSinglePoptRun=simSet.splitSinglePoptRun_;
    uint64_t startingTF0val = simSet.startingTF0val_, endTF0val=simSet.endTF0val_;
    
    

    long double popMeanFitness=zero;
    long double Popt=simSet.Popt_, omega=simSet.omega_;
    SimplestRegPathIndividual indivAABB(false), indivAABb(false), indivAAbb(false);//recombinants
    SimplestRegPathIndividual indivAaBB(false), indivAaBb(false), indivAabb(false);//
    SimplestRegPathIndividual indivaaBB(false), indivaaBb(false), indivaabb(false);
    long double wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,wBarMax;
    long double phAABB,phAABb,phAAbb,phAaBB,phAaBb,phAabb,phaaBB,phaaBb,phaabb;
    wAABB=wAABb=wAAbb=wAaBB=wAaBb=wAabb=waaBB=waaBb=waabb=wBarMax=popMeanFitness=zero;//fitnesses
    phAABB=phAABb=phAAbb=phAaBB=phAaBb=phAabb=phaaBB=phaaBb=phaabb=zero;//phenotypes
    std::vector<long double> phat(10,-one), qhat(10,-one),popMeanPhenotypes(10,-one);
    std::vector<bool> pNeutral(4,false), qNeutral(4,false);
    for(int i=0;i<4;++i){pNeutral[i]=qNeutral[i]=false;}
    wAABB=wAABb=wAAbb=wAaBB=wAaBb=wAabb=waaBB=waaBb=waabb=-one;
    phAABB=phAABb=phAAbb=phAaBB=phAaBb=phAabb=phaaBB=phaaBb=phaabb=-one;
    focalIndiv.CalculatePhenotype(simSet);
    focalIndiv.CalculateFitness(simSet.Popt_,simSet.omega_);
    std::string focalGtypeBitstringStr=focalIndiv.gtypeString(simSet.bitstringLen_);
    std::string focalGtypeMismatchStr=focalIndiv.mismatchStringMathematicaFormat();
    //create recombinant genotypes
    int numMaxima=0;
    if(focalIndiv.IsTFheterozygote() && focalIndiv.IsCisHeterozygote()){//maximize for p & q
            indivAaBb=focalIndiv;
            
            indivAABB.TFdosage_[0]=indivAABB.TFdosage_[1]=gtypeSet.dosageVal0_;
            indivAABB.TFproduct_[0]=indivAABB.TFproduct_[1]=gtypeSet.tfVal0_;
            indivAABB.cis_[0]=indivAABB.cis_[1]=gtypeSet.cisVal0_;
        
            indivAABb.TFdosage_[0]=indivAABb.TFdosage_[1]=gtypeSet.dosageVal0_;
            indivAABb.TFproduct_[0]=indivAABb.TFproduct_[1]=gtypeSet.tfVal0_;
            indivAABb.cis_[0]=gtypeSet.cisVal0_; indivAABb.cis_[1]=gtypeSet.cisVal1_;
        
            indivAAbb.TFdosage_[0]=indivAAbb.TFdosage_[1]=gtypeSet.dosageVal0_;
            indivAAbb.TFproduct_[0]=indivAAbb.TFproduct_[1]=gtypeSet.tfVal0_;
            indivAAbb.cis_[0]=indivAAbb.cis_[1]=gtypeSet.cisVal1_;
        
            indivAaBB.TFdosage_[0]=gtypeSet.dosageVal0_; indivAaBB.TFdosage_[1]=gtypeSet.dosageVal1_;
            indivAaBB.TFproduct_[0]=gtypeSet.tfVal0_; indivAaBB.TFproduct_[1]=gtypeSet.tfVal1_;
            indivAaBB.cis_[0]=indivAaBB.cis_[1]=gtypeSet.cisVal0_;
        
            indivAabb.TFdosage_[0]=gtypeSet.dosageVal0_; indivAabb.TFdosage_[1]=gtypeSet.dosageVal1_;
            indivAabb.TFproduct_[0]=gtypeSet.tfVal0_; indivAabb.TFproduct_[1]=gtypeSet.tfVal1_;
            indivAabb.cis_[0]=indivAabb.cis_[1]=gtypeSet.cisVal1_;
        
            indivaaBB.TFdosage_[0]=indivaaBB.TFdosage_[1]=gtypeSet.dosageVal1_;
            indivaaBB.TFproduct_[0]=indivaaBB.TFproduct_[1]=gtypeSet.tfVal1_;
            indivaaBB.cis_[0]=indivaaBB.cis_[1]=gtypeSet.cisVal0_;
        
            indivaaBb.TFdosage_[0]=indivaaBb.TFdosage_[1]=gtypeSet.dosageVal1_;
            indivaaBb.TFproduct_[0]=indivaaBb.TFproduct_[1]=gtypeSet.tfVal1_;
            indivaaBb.cis_[0]=gtypeSet.cisVal0_; indivaaBb.cis_[1]=gtypeSet.cisVal1_;
        
            indivaabb.TFdosage_[0]=indivaabb.TFdosage_[1]=gtypeSet.dosageVal1_;
            indivaabb.TFproduct_[0]=indivaabb.TFproduct_[1]=gtypeSet.tfVal1_;
            indivaabb.cis_[0]=indivaabb.cis_[1]=gtypeSet.cisVal1_;


            indivAABB.Reset();indivAABb.Reset();indivAAbb.Reset();indivAaBB.Reset();
            indivAabb.Reset();indivaaBB.Reset();indivaaBb.Reset();indivaabb.Reset();
            std::string focal=focalIndiv.mismatchStringMathematicaFormat();
            std::string mAABB=indivAABB.mismatchStringMathematicaFormat();
            std::string mAABb=indivAABb.mismatchStringMathematicaFormat();
            std::string mAAbb=indivAAbb.mismatchStringMathematicaFormat();
            std::string mAaBB=indivAaBB.mismatchStringMathematicaFormat();
            std::string mAaBb=focalIndiv.mismatchStringMathematicaFormat();
            std::string mAabb=indivAabb.mismatchStringMathematicaFormat();
            std::string maaBB=indivaaBB.mismatchStringMathematicaFormat();
            std::string maaBb=indivaaBb.mismatchStringMathematicaFormat();
            std::string maabb=indivaabb.mismatchStringMathematicaFormat();
            indivAABB.CalculatePhenotype(simSet);
            indivAABb.CalculatePhenotype(simSet);
            indivAAbb.CalculatePhenotype(simSet);
            indivAaBB.CalculatePhenotype(simSet);
            indivAabb.CalculatePhenotype(simSet);
            indivaaBB.CalculatePhenotype(simSet);
            indivaaBb.CalculatePhenotype(simSet);
            indivaabb.CalculatePhenotype(simSet);

            phAABB=indivAABB.phenotype();
            phAABb=indivAABb.phenotype();
            phAAbb=indivAAbb.phenotype();
            phAaBB=indivAaBB.phenotype();
            phAaBb=focalIndiv.phenotype();
            phAabb=indivAabb.phenotype();
            phaaBB=indivaaBB.phenotype();
            phaaBb=indivaaBb.phenotype();
            phaabb=indivaabb.phenotype();

            wAABB=indivAABB.CalculateFitness(Popt,omega);
            wAABb=indivAABb.CalculateFitness(Popt,omega);
            wAAbb=indivAAbb.CalculateFitness(Popt,omega);
            wAaBB=indivAaBB.CalculateFitness(Popt,omega);
            wAaBb=focalIndiv.fitness();
            wAabb=indivAabb.CalculateFitness(Popt,omega);
            waaBB=indivaaBB.CalculateFitness(Popt,omega);
            waaBb=indivaaBb.CalculateFitness(Popt,omega);
            waabb=indivaabb.CalculateFitness(Popt,omega);
        
        //maximize for p & q
        MaximizePopMeanFitnessPandQv2(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,popMeanFitness,phat,qhat,pNeutral,qNeutral,numMaxima,focalGtypeBitstringStr);
        collectMeanPhenotypesPQ(phAABB,phAABb,phAAbb,phAaBB,phAaBb,phAabb,phaaBB,phaaBb,phaabb,phat,qhat,popMeanPhenotypes,numMaxima);
            }
        else if(focalIndiv.IsTFheterozygote()){//maximize for p
            indivAaBB=focalIndiv;
            
            indivAABB.TFdosage_[0]=indivAABB.TFdosage_[1]=gtypeSet.dosageVal0_;
            indivAABB.TFproduct_[0]=indivAABB.TFproduct_[1]=gtypeSet.tfVal0_;
            indivAABB.cis_[0]=indivAABB.cis_[1]=gtypeSet.cisVal0_;

            indivaaBB.TFdosage_[0]=indivaaBB.TFdosage_[1]=gtypeSet.dosageVal1_;
            indivaaBB.TFproduct_[0]=indivaaBB.TFproduct_[1]=gtypeSet.tfVal1_;
            indivaaBB.cis_[0]=indivaaBB.cis_[1]=gtypeSet.cisVal0_;
            
            indivAABB.Reset();indivaaBB.Reset();
            indivAABB.CalculatePhenotype(simSet);
            indivaaBB.CalculatePhenotype(simSet);

            phAABB=indivAABB.phenotype();
            phAaBB=focalIndiv.phenotype();
            phaaBB=indivaaBB.phenotype();

            wAABB=indivAABB.CalculateFitness(Popt,omega);
            wAaBB=focalIndiv.fitness();
            waaBB=indivaaBB.CalculateFitness(Popt,omega);
            MaximizePopMeanFitnessP(wAABB,wAaBB,waaBB,popMeanFitness,phat,pNeutral,numMaxima);
            for(int m=0;m<numMaxima;++m){qhat[m]=one;qNeutral[m]=false;}
            collectMeanPhenotypesP(phAABB,phAaBB,phaaBB,phat,popMeanPhenotypes,numMaxima);
            }
        else if(focalIndiv.IsCisHeterozygote()){//maximize for q
            indivAABb = focalIndiv;
            
            indivAABB.TFdosage_[0]=indivAABB.TFdosage_[1]=gtypeSet.dosageVal0_;
            indivAABB.TFproduct_[0]=indivAABB.TFproduct_[1]=gtypeSet.tfVal0_;
            indivAABB.cis_[0]=indivAABB.cis_[1]=gtypeSet.cisVal0_;

            indivAAbb.TFdosage_[0]=indivAAbb.TFdosage_[1]=gtypeSet.dosageVal0_;
            indivAAbb.TFproduct_[0]=indivAAbb.TFproduct_[1]=gtypeSet.tfVal0_;
            indivAAbb.cis_[0]=indivAAbb.cis_[1]=gtypeSet.cisVal1_;

            indivAABB.Reset();indivAAbb.Reset();
            indivAABB.CalculatePhenotype(simSet);
            indivAAbb.CalculatePhenotype(simSet);
            phAABB=indivAABB.phenotype();
            phAABb=focalIndiv.phenotype();
            phAAbb=indivAAbb.phenotype();

            wAABB=indivAABB.CalculateFitness(Popt,omega);
            wAABb=focalIndiv.fitness();
            wAAbb=indivAAbb.CalculateFitness(Popt,omega);
            MaximizePopMeanFitnessQ(wAABB,wAABb,wAAbb,popMeanFitness,qhat,qNeutral,numMaxima);
            for(int m=0;m<numMaxima;++m){phat[m]=one;pNeutral[m]=false;}
            collectMeanPhenotypesQ(phAABB,phAABb,phAAbb,qhat,popMeanPhenotypes,numMaxima);
            }
        else{//double homozygote AABB
            indivAABB=focalIndiv;
            phAABB=focalIndiv.phenotype();
            wAABB=focalIndiv.fitness();
            phat[0]=qhat[0]=one;
            pNeutral[0]=qNeutral[0]=false;
            numMaxima=1;
            popMeanFitness=focalIndiv.fitness();
            popMeanPhenotypes[0]=focalIndiv.phenotype();
            }

    popMeanFitness=MIN(MAX(ROUND(popMeanFitness,decimalDigitsToRound),zero),one);
    for(int m=0;m<numMaxima;++m){
        popMeanPhenotypes[m]=MIN(MAX(ROUND(popMeanPhenotypes[m],decimalDigitsToRound),zero),one);}

    
    
    
    if(popMeanFitness==maxPopMeanFitness){

        for(int m=0;m<numMaxima;++m){
            if(phat[m]==-one && qhat[m]==-one) break;//reached the end
            std::string ht("___"), phatstr, qhatstr;//, ABgtype;
            std::stringstream phatstrSS,qhatstrSS;
            phatstrSS<<phat[m]; qhatstrSS<<qhat[m];
//                                    phatstr=std::to_string(phat[m]); qhatstr=std::to_string(qhat[m]);
            phatstr=phatstrSS.str(); qhatstr=phatstrSS.str();
            SimplestRegPathIndividual solutionIndiv;
            BitstringGenotypeData referenceIndivData, solutionIndivData;
            if(phat[m]==one){//the A TF allele is fixed
                    if(qhat[m]==one){//the B cis allele is fixed
                            solutionIndiv=indivAABB;
                            }
                        else if(qhat[m]==zero){//the b cis allele is fixed
                            solutionIndiv=indivAAbb;
                            }
                        else if(qNeutral[m]){
                            qhatstr="neutral";
                            ht[2]='n';
                            solutionIndiv=indivAABb;
                            }
                        else{//cis locus is heterozygous
                            ht[2]='c';
                            solutionIndiv=indivAABb;
                            }
                    }
                else if(phat[m]==zero){//the a TF allele is fixed
                    if(qhat[m]==one){//the B cis allele is fixed
                            solutionIndiv=indivaaBB;
                            }
                        else if(qhat[m]==zero){//the b cis allele is fixed
                            solutionIndiv=indivaabb;
                            }
                        else if(qNeutral[m]){
                            qhatstr="neutral";
                            ht[2]='n';
                            solutionIndiv=indivaaBb;
                            }
                        else{//cis locus is heterozygous
                            ht[2]='c';
                            solutionIndiv=indivaaBb;
                            }
                    }
                else if(pNeutral[m]){
                     phatstr="neutral";
                    //figure out which piece of the TF is heterozygous, if either
                    if(gtypeSet.dosageVal0_ != gtypeSet.dosageVal1_){
                        ht[0]='n';}//dosages differ
                    if(gtypeSet.tfVal0_ != gtypeSet.tfVal1_){
                        ht[1]='n';}//product sites differ
                    if(qhat[m]==one){//the B cis allele is fixed
                            solutionIndiv=indivAaBB;
                            }
                        else if(qhat[m]==zero){//the b cis allele is fixed
                            solutionIndiv=indivAabb;
                            }
                        else if(qNeutral[m]){
                            qhatstr="neutral";
                            ht[2]='n';
                            solutionIndiv=focalIndiv;//AaBb
                            }
                        else{//cis locus is heterozygous
                            ht[2]='c';
                            solutionIndiv=focalIndiv;//AaBb
                            }
                    }//neutral TF
                else{//the TF is polymorphic & non-neutral
                    //figure out which pieces of the TF are heterozygous/neutral
                    if(gtypeSet.dosageVal0_ != gtypeSet.dosageVal1_){
                        ht[0]='d';}//dosages differ and aren't neutral
                    if(qhat[m]==one){//the B cis allele is fixed
                            solutionIndiv=indivAaBB;
                            }
                        else if(qhat[m]==zero){//the b cis allele is fixed
                            solutionIndiv=indivAabb;
                            }
                        else if(qNeutral[m]){
                            ht[2]='n';
                            qhatstr="neutral";
                            solutionIndiv=focalIndiv;//AaBb
                            }
                        else{//cis locus is heterozygous
                            ht[2]='c';
                            solutionIndiv=focalIndiv;//AaBb
                            }
                    if(gtypeSet.tfVal0_ != gtypeSet.tfVal1_){//product sites differ
                        solutionIndiv.SetMismatchesUsingBitstrings();
                        if(solutionIndiv.IsTFprodHeterozygote(true)){//mismatch heterozygote
                                ht[1]='p';}
                            else{
                                ht[1]='n';}
                            }
                    }//polymorphic TF

            referenceIndivData.CollectData(focalIndiv,simSet.bitstringLen_);
            solutionIndiv.CalculatePhenotype(simSet);//sets mismatch pattern
            //see if this solution has been found via another reference sequence
            long double pMostCommon=-one,qMostCommon=-one;
            if(!pNeutral[m]){pMostCommon=MAX(phat[m],one-phat[m]);}
            if(!qNeutral[m]){qMostCommon=MAX(qhat[m],one-qhat[m]);}
            solutionIndiv.SetMismatchesUsingBitstrings();
            std::string solGtype=solutionIndiv.gtypeString(simSet.bitstringLen_);
            std::string refGtype=focalIndiv.gtypeString(simSet.bitstringLen_);
            std::string mhp = solutionIndiv.mismatchHetType();
            std::string mp = solutionIndiv.mismatchStringMathematicaFormat();
            int hc = codeToInt(ht), mhc = codeToInt(mhp);
            FitnessMaximumSolutionSet thisSolutionSummary(Popt,omega,simSet.NtfsatPerAllele_,simSet.bitstringLen_,popMeanFitness,
                                                          popMeanPhenotypes[m],pMostCommon,
                                                          qMostCommon,pNeutral[m],qNeutral[m],
                                                          hc,ht,mhc,mhp,mp,solGtype,refGtype,
                                                          splitSinglePoptRun,startingTF0val,endTF0val);
            fmss=thisSolutionSummary;
            solutionEqualsMaxOrBetter=true;
            }//m
        }//popMeanFitness==maxPopMeanFitness

    if(popMeanFitness>maxPopMeanFitness){
//        newBest=true;
        for(int m=0;m<numMaxima;++m){
            SimplestRegPathIndividual solutionIndiv;
            BitstringGenotypeData referenceIndivData, solutionIndivData;

            std::string ht("___"), phatstr, qhatstr;
            phatstr=std::to_string(phat[m]); qhatstr=std::to_string(qhat[m]);
            if(phat[m]==one){//the A TF allele is fixed
                    if(qhat[m]==one){//the B cis allele is fixed
                            solutionIndiv=indivAABB;
                            }
                        else if(qhat[m]==zero){//the b cis allele is fixed
                            solutionIndiv=indivAAbb;
                            }
                        else if(qNeutral[m]){
                            qhatstr="neutral";
                            ht[2]='n';
                            solutionIndiv=indivAABb;
                            }
                        else{//cis locus is heterozygous
                            ht[2]='c';
                            solutionIndiv=indivAABb;
                            }
                    }
                else if(phat[m]==zero){//the a TF allele is fixed
                    if(qhat[m]==one){//the B cis allele is fixed
                            solutionIndiv=indivaaBB;
                            }
                        else if(qhat[m]==zero){//the b cis allele is fixed
                            solutionIndiv=indivaabb;
                            }
                        else if(qNeutral[m]){
                            qhatstr="neutral";
                            ht[2]='n';
                            solutionIndiv=indivaaBb;
                            }
                        else{//cis locus is heterozygous
                            ht[2]='c';
                            solutionIndiv=indivaaBb;
                            }
                    }
                else if(pNeutral[m]){
                    phatstr="neutral";
                    //figure out which piece of the TF is heterozygous, if either
                    if(gtypeSet.dosageVal0_ != gtypeSet.dosageVal1_){
                        ht[0]='n';}//dosages differ; neutrality is here
                    if(gtypeSet.tfVal0_ != gtypeSet.tfVal1_){
                        ht[1]='n';}//product sites differ; neutrality is here
                    if(qhat[m]==one){//the B cis allele is fixed
                            solutionIndiv=indivAaBB;
                            }
                        else if(qhat[m]==zero){//the b cis allele is fixed
                            solutionIndiv=indivAabb;
                            }
                        else if(qNeutral[m]){//shouldn't get here
                            qhatstr="neutral";
                            ht[2]='n';
                            solutionIndiv=focalIndiv;
                            }
                        else{//cis locus is heterozygous
                            ht[2]='c';
                            solutionIndiv=focalIndiv;
                            }
                    }//neutral TF
                else{//the TF is polymorphic & non-neutral
                    //figure out which pieces of the TF are heterozygous/neutral
                    if(gtypeSet.dosageVal0_ != gtypeSet.dosageVal1_){
                        ht[0]='d';}//dosages differ and aren't neutral
                    if(qhat[m]==one){//the B cis allele is fixed
                            solutionIndiv=indivAaBB;
                            }
                        else if(qhat[m]==zero){//the b cis allele is fixed
                            solutionIndiv=indivAabb;
                            }
                        else if(qNeutral[m]){
                            ht[2]='n';
                            qhatstr="neutral";
                            solutionIndiv=focalIndiv;//AaBb
                            }
                        else{//cis locus is heterozygous
                            ht[2]='c';
                            solutionIndiv=focalIndiv;//AaBb
                            }
                    if(gtypeSet.tfVal0_ != gtypeSet.tfVal1_){//product sites differ
                        solutionIndiv.SetMismatchesUsingBitstrings();
                        if(solutionIndiv.IsTFprodHeterozygote(true)){//mismatch heterozygote
                                ht[1]='p';}
                            else{
                                ht[1]='n';}
                            }
                    }//polymorphic TF
            maxPopMeanFitness=popMeanFitness;
            referenceIndivData.CollectData(focalIndiv,simSet.bitstringLen_);
            solutionIndiv.CalculatePhenotype(simSet);//sets mismatch pattern data
            //see if this solution has been found via another reference sequence
            long double pMostCommon=-one,qMostCommon=-one;
            if(!pNeutral[m]){pMostCommon=MAX(phat[m],one-phat[m]);}
            if(!qNeutral[m]){qMostCommon=MAX(qhat[m],one-qhat[m]);}
            solutionIndiv.SetMismatchesUsingBitstrings();
            std::string solGtype=solutionIndiv.gtypeString(simSet.bitstringLen_);
            std::string refGtype=focalIndiv.gtypeString(simSet.bitstringLen_);
            std::string mhp = solutionIndiv.mismatchHetType();
            std::string mp = solutionIndiv.mismatchStringMathematicaFormat();
            int hc = codeToInt(ht), mhc = codeToInt(mhp);
            FitnessMaximumSolutionSet thisSolutionSummary(Popt,omega,simSet.NtfsatPerAllele_,simSet.bitstringLen_,popMeanFitness,
                                                          popMeanPhenotypes[m],pMostCommon,
                                                          qMostCommon,pNeutral[m],qNeutral[m],
                                                          hc,ht,mhc,mhp,mp,solGtype,refGtype,
                                                          splitSinglePoptRun,startingTF0val,endTF0val);
            fmss=thisSolutionSummary;
            solutionEqualsMaxOrBetter=true;
            }//m
        }//popMeanFitness>maxPopMeanFitness
    
    }//MaximizeUsingBitstringsOneReferenceGtype




void MaximizeUsingBitstringsDosageOnly(int bitstringLen, long double NtfsatPerAllele, long double deltaG1dosage,
                long double deltaG1, long double minExpression, long double maxExpression,
                long double Popt, long double omega,
                FitnessMaximaBitstringSolutions& wBarMaxAllSolutions, std::ostream& outputfileAllSolutions, bool saveAllSolutions,
                FitnessMaximaSolutionSets& summariesOfSolutions, std::ostream& outputfileSolutionSummaries){
    //this version considers only the dosage site
    //each allele 2 variant runs in its own thread
    coutLock.lock(); std::cout<<"maximizing for Popt="<<Popt<<" and omega="<<omega<<std::endl; coutLock.unlock();
    uint64_t maxBitstringVal = uint64_t(pow(2,bitstringLen));
    
    simulationSettings simSet(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression,Popt,omega);
    genotypeSettings gtypeSet;
    SimplestRegPathIndividual focalIndiv(false);//indivAaBb
    gtypeSet.tfVal0_=gtypeSet.tfVal1_=gtypeSet.cisVal0_=gtypeSet.cisVal1_=0;
    focalIndiv.SetGenotype(1,0,0); focalIndiv.SetGenotype(1,1,0);
    focalIndiv.SetGenotype(2,0,0); focalIndiv.SetGenotype(2,1,0);
    long double maxPopMeanFitness=-one;
    for(uint64_t dosageVal0=0;dosageVal0<maxBitstringVal;++dosageVal0){//1st TF allele dosage
        gtypeSet.dosageVal0_=dosageVal0;
        focalIndiv.SetGenotype(0,0,dosageVal0);
        uint64_t dosageVal1counter=0;
        std::thread refThread;
        while(dosageVal1counter<=dosageVal0){
            int availableThreads = MAX((unsigned int)1,refThread.hardware_concurrency());
            int threadsToUse = MIN(availableThreads,(int)(dosageVal0-dosageVal1counter)+1);//+1 for when disageVal0==dosageVal1
            threadsToUse = MIN(threadsToUse,16);
            if(threadsToUse<=1){//just call the function without using a thread
                    gtypeSet.dosageVal1_=dosageVal1counter;
                    focalIndiv.SetGenotype(0,1,dosageVal1counter);
                    int solutionEqualsMaxOrBetter=false;
                    FitnessMaximumSolutionSet newSolutionSummary;
                    MaximizeUsingBitstringsOneReferenceGtype(&focalIndiv,&simSet,&gtypeSet,&newSolutionSummary,
                        maxPopMeanFitness,&solutionEqualsMaxOrBetter);
                    if(solutionEqualsMaxOrBetter){
                        maxPopMeanFitness=newSolutionSummary.wBarMax();
                        summariesOfSolutions.AddSolution(newSolutionSummary);}
                    dosageVal1counter++;
                    }
                else{//make & use threads
                    std::vector<std::thread> threadList;
                    std::vector<int> solutionSameOrBetterTF;
                    std::vector<FitnessMaximumSolutionSet> newSolutionSummaries;
                    std::vector<SimplestRegPathIndividual> focalIndivsToTest;
                    std::vector<genotypeSettings> gtypeSettingsList;
                    std::vector<simulationSettings> simSettingsList;
                    for(int t=0;t<threadsToUse;++t){//initialize
                        solutionSameOrBetterTF.push_back(false);
                        FitnessMaximumSolutionSet fmss;
                        newSolutionSummaries.push_back(fmss);
                        gtypeSet.dosageVal1_=dosageVal1counter;
                        gtypeSettingsList.push_back(gtypeSet);
                        simSettingsList.push_back(simSet);
                        focalIndiv.SetGenotype(0,1,dosageVal1counter);
                        focalIndivsToTest.push_back(focalIndiv);
                        dosageVal1counter++;
                        }//initialize
                    for(int t=0;t<threadsToUse;++t){//make and run threads
                        SimplestRegPathIndividual *fIndiv = &(focalIndivsToTest[t]);
                        FitnessMaximumSolutionSet *fmss = &(newSolutionSummaries[t]);
                        genotypeSettings *gset = &(gtypeSettingsList[t]);
                        simulationSettings *simset = &(simSettingsList[t]);
                        int *sameOrBetter = &(solutionSameOrBetterTF[t]);
                        std::thread th(
                            MaximizeUsingBitstringsOneReferenceGtype,fIndiv,simset,gset,fmss,
                                            maxPopMeanFitness,sameOrBetter);
                        threadList.push_back(move(th));
                        }//t
                    for(int t=0;t<threadList.size();++t){//synchronize threads
                        threadList[t].join();}
                    for(int t=0;t<threadList.size();++t){//collect data
                        if(solutionSameOrBetterTF[t]){
                            summariesOfSolutions.AddSolution(newSolutionSummaries[t]);
                            maxPopMeanFitness=MAX(maxPopMeanFitness,newSolutionSummaries[t].wBarMax());
                            }//solutionSameOrBetterTF
                        }//t
                    }//make & use threads
            }//while dosageVal1counter
        }//dosageVal0
    coutLock.lock(); std::cout<<std::endl<<"**************************"<<std::endl; coutLock.unlock();
    summariesOfSolutions.PrintDataByPopt(outputfileSolutionSummaries, Popt);
    }//MaximizeUsingBitstringsDosageOnly


void MaximizeUsingBitstringsTFproductOnly(int bitstringLen, long double NtfsatPerAllele, long double deltaG1dosage,
                long double deltaG1, long double minExpression, long double maxExpression,
                long double Popt, long double omega,
                FitnessMaximaBitstringSolutions& wBarMaxAllSolutions, std::ostream& outputfileAllSolutions, bool saveAllSolutions,
                FitnessMaximaSolutionSets& summariesOfSolutions, std::ostream& outputfileSolutionSummaries){
    //this version considers only the dosage site
    //each allele 2 variant runs in its own thread
    coutLock.lock(); std::cout<<"maximizing for Popt="<<Popt<<" and omega="<<omega<<std::endl; coutLock.unlock();
    uint64_t maxBitstringVal = uint64_t(pow(2,bitstringLen));
    
    simulationSettings simSet(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression,Popt,omega);
    genotypeSettings gtypeSet;
    SimplestRegPathIndividual focalIndiv(false);//indivAaBb
    gtypeSet.dosageVal0_=gtypeSet.dosageVal1_=gtypeSet.cisVal0_=gtypeSet.cisVal1_=0;
    focalIndiv.SetGenotype(0,0,0); focalIndiv.SetGenotype(0,1,0);
    focalIndiv.SetGenotype(2,0,0); focalIndiv.SetGenotype(2,1,0);
    long double maxPopMeanFitness=-one;
    for(uint64_t tfVal0=0;tfVal0<maxBitstringVal;++tfVal0){//1st TF allele product
        gtypeSet.tfVal0_=tfVal0;
        focalIndiv.SetGenotype(1,0,tfVal0);
        uint64_t tfVal1counter=0;
        std::thread refThread;
        while(tfVal1counter<=tfVal0){
            int availableThreads = MAX((unsigned int)1,refThread.hardware_concurrency());
            int threadsToUse = MIN(availableThreads,(int)(tfVal0-tfVal1counter)+1);//+1 for when disageVal0==dosageVal1
            threadsToUse = MIN(threadsToUse,16);
            if(threadsToUse<=1){//just call the function without using a thread
                    gtypeSet.tfVal1_=tfVal1counter;
                    focalIndiv.SetGenotype(1,1,tfVal1counter);
                    int solutionEqualsMaxOrBetter=false;
                    FitnessMaximumSolutionSet newSolutionSummary;
                    MaximizeUsingBitstringsOneReferenceGtype(&focalIndiv,&simSet,&gtypeSet,&newSolutionSummary,
                        maxPopMeanFitness,&solutionEqualsMaxOrBetter);
                    if(solutionEqualsMaxOrBetter){
                        maxPopMeanFitness=newSolutionSummary.wBarMax();
                        summariesOfSolutions.AddSolution(newSolutionSummary);}
                    tfVal1counter++;
                    }
                else{//make & use threads
                    std::vector<std::thread> threadList;
                    std::vector<int> solutionSameOrBetterTF;
                    std::vector<FitnessMaximumSolutionSet> newSolutionSummaries;
                    std::vector<SimplestRegPathIndividual> focalIndivsToTest;
                    std::vector<genotypeSettings> gtypeSettingsList;
                    std::vector<simulationSettings> simSettingsList;
                    for(int t=0;t<threadsToUse;++t){//initialize
                        solutionSameOrBetterTF.push_back(false);
                        FitnessMaximumSolutionSet fmss;
                        newSolutionSummaries.push_back(fmss);
                        gtypeSet.tfVal1_=tfVal1counter;
                        gtypeSettingsList.push_back(gtypeSet);
                        simSettingsList.push_back(simSet);
                        focalIndiv.SetGenotype(1,1,tfVal1counter);
                        focalIndivsToTest.push_back(focalIndiv);
                        tfVal1counter++;
                        }//initialize
                    for(int t=0;t<threadsToUse;++t){//make and run threads
                        SimplestRegPathIndividual *fIndiv = &(focalIndivsToTest[t]);
                        FitnessMaximumSolutionSet *fmss = &(newSolutionSummaries[t]);
                        genotypeSettings *gset = &(gtypeSettingsList[t]);
                        simulationSettings *simset = &(simSettingsList[t]);
                        int *sameOrBetter = &(solutionSameOrBetterTF[t]);
                        std::thread th(
                            MaximizeUsingBitstringsOneReferenceGtype,fIndiv,simset,gset,fmss,
                                            maxPopMeanFitness,sameOrBetter);
                        threadList.push_back(move(th));
                        }//t
                    for(int t=0;t<threadList.size();++t){//synchronize threads
                        threadList[t].join();}
                    for(int t=0;t<threadList.size();++t){//collect data
                        if(solutionSameOrBetterTF[t]){
                            summariesOfSolutions.AddSolution(newSolutionSummaries[t]);
                            maxPopMeanFitness=MAX(maxPopMeanFitness,newSolutionSummaries[t].wBarMax());
                            }//solutionSameOrBetterTF
                        }//t
                    }//make & use threads
            }//while tfVal1counter
        }//tfVal0
    coutLock.lock(); std::cout<<std::endl<<"**************************"<<std::endl; coutLock.unlock();
    summariesOfSolutions.PrintDataByPopt(outputfileSolutionSummaries, Popt);
    }//MaximizeUsingBitstringsTFproductOnly



void MaximizeUsingBitstringsCisOnly(int bitstringLen, long double NtfsatPerAllele, long double deltaG1dosage,
                long double deltaG1, long double minExpression, long double maxExpression,
                long double Popt, long double omega,
                FitnessMaximaBitstringSolutions& wBarMaxAllSolutions, std::ostream& outputfileAllSolutions, bool saveAllSolutions,
                FitnessMaximaSolutionSets& summariesOfSolutions, std::ostream& outputfileSolutionSummaries){
    //this version considers only the cis site
    //each allele 2 variant runs in its own thread
    coutLock.lock(); std::cout<<"maximizing for Popt="<<Popt<<" and omega="<<omega<<std::endl; coutLock.unlock();
    uint64_t maxBitstringVal = uint64_t(pow(2,bitstringLen));
    
    simulationSettings simSet(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression,Popt,omega);
    genotypeSettings gtypeSet;
    SimplestRegPathIndividual focalIndiv(false);//indivAaBb
    gtypeSet.dosageVal0_=gtypeSet.dosageVal1_=gtypeSet.tfVal0_=gtypeSet.tfVal1_=0;
    focalIndiv.SetGenotype(0,0,0); focalIndiv.SetGenotype(0,1,0);
    focalIndiv.SetGenotype(1,0,0); focalIndiv.SetGenotype(1,1,0);
    long double maxPopMeanFitness=-one;
    for(uint64_t cisVal0=0;cisVal0<maxBitstringVal;++cisVal0){//1st TF allele product
        gtypeSet.cisVal0_=cisVal0;
        focalIndiv.SetGenotype(2,0,cisVal0);
        uint64_t cisVal1counter=0;
        std::thread refThread;
        while(cisVal1counter<=cisVal0){
            int availableThreads = MAX((unsigned int)1,refThread.hardware_concurrency());
            int threadsToUse = MIN(availableThreads,(int)(cisVal0-cisVal1counter)+1);//+1 for when disageVal0==dosageVal1
            threadsToUse = MIN(threadsToUse,16);
            if(threadsToUse<=1){//just call the function without using a thread
                    gtypeSet.tfVal1_=cisVal1counter;
                    focalIndiv.SetGenotype(2,1,cisVal1counter);
                    int solutionEqualsMaxOrBetter=false;
                    FitnessMaximumSolutionSet newSolutionSummary;
                    MaximizeUsingBitstringsOneReferenceGtype(&focalIndiv,&simSet,&gtypeSet,&newSolutionSummary,
                        maxPopMeanFitness,&solutionEqualsMaxOrBetter);
                    if(solutionEqualsMaxOrBetter){
                        maxPopMeanFitness=newSolutionSummary.wBarMax();
                        summariesOfSolutions.AddSolution(newSolutionSummary);}
                    cisVal1counter++;
                    }
                else{//make & use threads
                    std::vector<std::thread> threadList;
                    std::vector<int> solutionSameOrBetterTF;
                    std::vector<FitnessMaximumSolutionSet> newSolutionSummaries;
                    std::vector<SimplestRegPathIndividual> focalIndivsToTest;
                    std::vector<genotypeSettings> gtypeSettingsList;
                    std::vector<simulationSettings> simSettingsList;
                    for(int t=0;t<threadsToUse;++t){//initialize
                        solutionSameOrBetterTF.push_back(false);
                        FitnessMaximumSolutionSet fmss;
                        newSolutionSummaries.push_back(fmss);
                        gtypeSet.cisVal1_=cisVal1counter;
                        gtypeSettingsList.push_back(gtypeSet);
                        simSettingsList.push_back(simSet);
                        focalIndiv.SetGenotype(2,1,cisVal1counter);
                        focalIndivsToTest.push_back(focalIndiv);
                        cisVal1counter++;
                        }//initialize
                    for(int t=0;t<threadsToUse;++t){//make and run threads
                        SimplestRegPathIndividual *fIndiv = &(focalIndivsToTest[t]);
                        FitnessMaximumSolutionSet *fmss = &(newSolutionSummaries[t]);
                        genotypeSettings *gset = &(gtypeSettingsList[t]);
                        simulationSettings *simset = &(simSettingsList[t]);
                        int *sameOrBetter = &(solutionSameOrBetterTF[t]);
                        std::thread th(
                            MaximizeUsingBitstringsOneReferenceGtype,fIndiv,simset,gset,fmss,
                                            maxPopMeanFitness,sameOrBetter);
                        threadList.push_back(move(th));
                        }//t
                    for(int t=0;t<threadList.size();++t){//synchronize threads
                        threadList[t].join();}
                    for(int t=0;t<threadList.size();++t){//collect data
                        if(solutionSameOrBetterTF[t]){
                            summariesOfSolutions.AddSolution(newSolutionSummaries[t]);
                            maxPopMeanFitness=MAX(maxPopMeanFitness,newSolutionSummaries[t].wBarMax());
                            }//solutionSameOrBetterTF
                        }//t
                    }//make & use threads
            }//while cisVal1counter
        }//cisVal0
    coutLock.lock(); std::cout<<std::endl<<"**************************"<<std::endl; coutLock.unlock();
    summariesOfSolutions.PrintDataByPopt(outputfileSolutionSummaries, Popt);
    }//MaximizeUsingBitstringsCisOnly





void MaximizeUsingBitstringsAllTF1Gtypes(SimplestRegPathIndividual* focalIndivP, simulationSettings* simSetP,
        genotypeSettings* gtypeSetP, FitnessMaximaSolutionSets *summariesOfSolutionsP, long double maxPopMeanFitness,
        int* solutionEqualsMaxOrBetterP){
    //for a given full TF genotype at allele 0 and TF dosage genotype at allele 1,
    //this checks all TF-structural genotypes for allele 2 to find the one with the best fitness
    //the cis locus is held constant so that longer bitstrings can be studied
    SimplestRegPathIndividual& focalIndiv = *focalIndivP;
    simulationSettings& simSet = *simSetP;
    genotypeSettings& gtypeSet = *gtypeSetP;
    FitnessMaximaSolutionSets& summariesOfSolutions = *summariesOfSolutionsP;
    int& solutionEqualsMaxOrBetter = *solutionEqualsMaxOrBetterP;
    FitnessMaximumSolutionSet newSolutionSummary;
    long double oldMaxPopMeanFitness = maxPopMeanFitness;
    uint64_t tfVal0 = gtypeSet.tfVal0_;
    for(uint64_t tfVal1=0;tfVal1<=tfVal0;++tfVal1){
        focalIndiv.SetGenotype(1,1,tfVal1);
        gtypeSet.tfVal1_=tfVal1;
        solutionEqualsMaxOrBetter=false;
        newSolutionSummary.Reset();
        MaximizeUsingBitstringsOneReferenceGtype(&focalIndiv,&simSet,&gtypeSet,&newSolutionSummary,
            maxPopMeanFitness,&solutionEqualsMaxOrBetter);
        if(solutionEqualsMaxOrBetter){
            maxPopMeanFitness=newSolutionSummary.wBarMax();
            summariesOfSolutions.AddSolution(newSolutionSummary);}
        }//tfVal1
    solutionEqualsMaxOrBetter=false;
    if(maxPopMeanFitness>=oldMaxPopMeanFitness){
        solutionEqualsMaxOrBetter=true;}
    }//MaximizeUsingBitstringsAllTF1Gtypes




void MaximizeUsingBitstringsTFOnly(int bitstringLen, long double NtfsatPerAllele, long double deltaG1dosage,
                long double deltaG1, long double minExpression, long double maxExpression,
                long double Popt, long double omega,
                FitnessMaximaBitstringSolutions& wBarMaxAllSolutions, std::ostream& outputfileAllSolutions, bool saveAllSolutions,
                FitnessMaximaSolutionSets& summariesOfSolutions, std::ostream& outputfileSolutionSummaries){
    //this version holds the cis locus constant, maximizing over variation in the TF genotype
    //each allele 2 variant runs in its own thread
    coutLock.lock(); std::cout<<"maximizing for Popt="<<Popt<<" and omega="<<omega<<std::endl; coutLock.unlock();
    uint64_t maxBitstringVal = uint64_t(pow(2,bitstringLen));
    
    simulationSettings simSet(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression,Popt,omega);
    genotypeSettings gtypeSet;
    SimplestRegPathIndividual focalIndiv(false);//indivAaBb
    gtypeSet.cisVal0_=gtypeSet.cisVal1_=0;
    focalIndiv.SetGenotype(2,0,0); focalIndiv.SetGenotype(2,1,0);
    long double maxPopMeanFitness=-one;
    for(uint64_t dosageVal0=0;dosageVal0<maxBitstringVal;++dosageVal0){//1st TF allele dosage
        gtypeSet.dosageVal0_=dosageVal0;
        focalIndiv.SetGenotype(0,0,dosageVal0);
        for(uint64_t tfVal0=0;tfVal0<maxBitstringVal;++tfVal0){//1st TF allele product
            gtypeSet.tfVal0_=tfVal0;
            focalIndiv.SetGenotype(1,0,tfVal0);
            uint64_t dosageVal1counter=0;
            std::thread refThread;
            while(dosageVal1counter<=dosageVal0){
                int availableThreads = MAX((unsigned int)1,refThread.hardware_concurrency());
                int threadsToUse = MIN(availableThreads,(int)(dosageVal0-dosageVal1counter)+1);//+1 for when disageVal0==dosageVal1
                threadsToUse = MIN(threadsToUse,16);
                if(threadsToUse<=1){//just call the function without using a thread
                        gtypeSet.dosageVal1_=dosageVal1counter;
                        focalIndiv.SetGenotype(0,1,dosageVal1counter);
                        int solutionEqualsMaxOrBetter=false;
                        FitnessMaximaSolutionSets newSolutionSummaries;
                        MaximizeUsingBitstringsAllTF1Gtypes(&focalIndiv,&simSet,&gtypeSet,&newSolutionSummaries,
                            maxPopMeanFitness,&solutionEqualsMaxOrBetter);
                        if(solutionEqualsMaxOrBetter){
                            maxPopMeanFitness=newSolutionSummaries.wBarMax(Popt);
                            summariesOfSolutions.ConcatenateSolutions(newSolutionSummaries);}
                        dosageVal1counter++;
                        }
                    else{//make & use threads
                        std::vector<std::thread> threadList;
                        std::vector<int> solutionSameOrBetterTF;
                        std::vector<FitnessMaximaSolutionSets> newSolutionSummarySets;
                        std::vector<SimplestRegPathIndividual> focalIndivsToTest;
                        std::vector<genotypeSettings> gtypeSettingsList;
                        std::vector<simulationSettings> simSettingsList;
                        for(int t=0;t<threadsToUse;++t){//initialize
                            solutionSameOrBetterTF.push_back(false);
                            FitnessMaximaSolutionSets fmss;
                            newSolutionSummarySets.push_back(fmss);
                            gtypeSet.dosageVal1_=dosageVal1counter;
                            gtypeSettingsList.push_back(gtypeSet);
                            simSettingsList.push_back(simSet);
                            focalIndiv.SetGenotype(1,1,dosageVal1counter);
                            focalIndivsToTest.push_back(focalIndiv);
                            dosageVal1counter++;
                            }//initialize
                        for(int t=0;t<threadsToUse;++t){//make and run threads
                            SimplestRegPathIndividual *fIndiv = &(focalIndivsToTest[t]);
                            FitnessMaximaSolutionSets *fmss = &(newSolutionSummarySets[t]);
                            genotypeSettings *gset = &(gtypeSettingsList[t]);
                            simulationSettings *simset = &(simSettingsList[t]);
                            int *sameOrBetter = &(solutionSameOrBetterTF[t]);
                            std::thread th(
                                MaximizeUsingBitstringsAllTF1Gtypes,fIndiv,simset,gset,fmss,
                                                maxPopMeanFitness,sameOrBetter);
                            threadList.push_back(move(th));
                            }//t
                        for(int t=0;t<threadList.size();++t){//synchronize threads
                            threadList[t].join();}
                        for(int t=0;t<threadList.size();++t){//collect data
                            if(solutionSameOrBetterTF[t]){
                                summariesOfSolutions.ConcatenateSolutions(newSolutionSummarySets[t]);
                                maxPopMeanFitness=MAX(maxPopMeanFitness,newSolutionSummarySets[t].wBarMax(Popt));
                                }//solutionSameOrBetterTF
                            }//t
                        }//make & use threads
                }//while dosageVal1counter
            }//tfVal0
        }//dosageVal0
    coutLock.lock(); std::cout<<std::endl<<"**************************"<<std::endl; coutLock.unlock();
    summariesOfSolutions.PrintDataByPopt(outputfileSolutionSummaries, Popt);
    }//MaximizeUsingBitstringsTFOnly



void MaximizeUsingBitstringsAllCisGtypes(SimplestRegPathIndividual* focalIndivP, simulationSettings* simSetP,
        genotypeSettings* gtypeSetP, FitnessMaximaSolutionSets *summariesOfSolutionsP, long double maxPopMeanFitness,
        int* solutionEqualsMaxOrBetterP){
    //for a given TF genotype, this checks all cis-locus genotypes to find the one with the best fitness
    //The parameters are set up so that it can be run in its own thread
    SimplestRegPathIndividual& focalIndiv = *focalIndivP;
    simulationSettings& simSet = *simSetP;
    genotypeSettings& gtypeSet = *gtypeSetP;
    FitnessMaximaSolutionSets& summariesOfSolutions = *summariesOfSolutionsP;
    int& solutionEqualsMaxOrBetter = *solutionEqualsMaxOrBetterP;
    uint64_t maxBitstringVal = uint64_t(pow(2,simSet.bitstringLen_));
    bool splitSinglePoptRun=summariesOfSolutions.splitSinglePoptRun_;
    uint64_t startingTF0val=summariesOfSolutions.startingTF0val_;
    uint64_t endTF0val=summariesOfSolutions.endTF0val_;
    FitnessMaximumSolutionSet newSolutionSummary(simSet.bitstringLen_,splitSinglePoptRun,startingTF0val,endTF0val);
    long double oldMaxPopMeanFitness = maxPopMeanFitness;
    for(uint64_t cisVal0=0;cisVal0<maxBitstringVal;++cisVal0){//1st cis allele promoter
        focalIndiv.SetGenotype(2,0,cisVal0);
        gtypeSet.cisVal0_=cisVal0;
        for(uint64_t cisVal1=0;cisVal1<=cisVal0;++cisVal1){//2nd cis allele promoter
            focalIndiv.SetGenotype(2,1,cisVal1);
            gtypeSet.cisVal1_=cisVal1;
            solutionEqualsMaxOrBetter=false;
            newSolutionSummary.Reset();
            MaximizeUsingBitstringsOneReferenceGtype(&focalIndiv,&simSet,&gtypeSet,&newSolutionSummary,
                maxPopMeanFitness,&solutionEqualsMaxOrBetter);
            if(solutionEqualsMaxOrBetter){
                maxPopMeanFitness=newSolutionSummary.wBarMax();
                summariesOfSolutions.AddSolution(newSolutionSummary);}
            }//cisVal1
        }//cisVal0
    if(maxPopMeanFitness>=oldMaxPopMeanFitness){
        solutionEqualsMaxOrBetter=true;}
    }//MaximizeUsingBitstringsAllCisGtypes




void MaximizeUsingBitstringsThreadable(int bitstringLen, long double NtfsatPerAllele, long double deltaG1dosage,
                long double deltaG1, long double minExpression, long double maxExpression,
                long double Popt, long double omega,
                FitnessMaximaBitstringSolutions& wBarMaxAllSolutions, std::ostream& outputfileAllSolutions, bool saveAllSolutions,
                FitnessMaximaSolutionSets& summariesOfSolutions, std::ostream& outputfileSolutionSummaries){
    coutLock.lock(); std::cout<<"maximizing for Popt="<<Popt<<" and omega="<<omega<<std::endl; coutLock.unlock();
    uint64_t maxBitstringVal = uint64_t(pow(2,bitstringLen));
    
    simulationSettings simSet(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression,Popt,omega);
    genotypeSettings gtypeSet;
    SimplestRegPathIndividual focalIndiv(false);//indivAaBb
    long double maxPopMeanFitness=-one;
//    std::vector<long double> phat(10,-one), qhat(10,-one),popMeanPhenotypes(10,-one);
//    std::vector<bool> pNeutral(4,false), qNeutral(4,false);
    
    for(uint64_t dosageVal0=0;dosageVal0<maxBitstringVal;++dosageVal0){//1st TF allele dosage
        gtypeSet.dosageVal0_=dosageVal0;
        focalIndiv.SetGenotype(0,0,dosageVal0);
        for(uint64_t tfVal0=0;tfVal0<maxBitstringVal;++tfVal0){//1st TF allele product
            gtypeSet.tfVal0_=tfVal0;
            focalIndiv.SetGenotype(1,0,tfVal0);
            for(uint64_t dosageVal1=0;dosageVal1<=dosageVal0;++dosageVal1){//2nd TF allele dosage
                gtypeSet.dosageVal1_=dosageVal1;
                focalIndiv.SetGenotype(0,1,dosageVal1);
                for(uint64_t tfVal1=0;tfVal1<=tfVal0;++tfVal1){//2nd TF allele product
                    gtypeSet.tfVal1_=tfVal1;
                    focalIndiv.SetGenotype(1,1,tfVal1);
                    //could put here a thread that maximized fitness over the whole cis parameter space for this TF pair
                    //so, if using threads, create a local FitnessMaximaSolutionSets object
                    
                    for(uint64_t cisVal0=0;cisVal0<maxBitstringVal;++cisVal0){//1st cis allele promoter
                        //add a while loop here for threading the cisVal1 calculations
                        //determine the number of available threads to set the loop bounds; this may vary dynamically
                            uint64_t cisVal1counter=0;
                            std::thread refThread;
                            while(cisVal1counter<=cisVal0){
                                int availableThreads = MAX((unsigned int)1,refThread.hardware_concurrency());
                                int threadsToUse = MIN(availableThreads,(int)(cisVal0-cisVal1counter)+1);//+1 for when cisVal0==cisVal1
                                threadsToUse = MIN(threadsToUse,16);
                                if(threadsToUse<=1){//just call the function without using a thread
                                        gtypeSet.cisVal1_=cisVal1counter;
                                        focalIndiv.SetGenotype(2,1,cisVal1counter);
                                        int solutionEqualsMaxOrBetter=false;
                                        FitnessMaximumSolutionSet newSolutionSummary;
                                        MaximizeUsingBitstringsOneReferenceGtype(&focalIndiv,&simSet,&gtypeSet,&newSolutionSummary,
                                            maxPopMeanFitness,&solutionEqualsMaxOrBetter);
                                        if(solutionEqualsMaxOrBetter){
                                            maxPopMeanFitness=newSolutionSummary.wBarMax();
                                            summariesOfSolutions.AddSolution(newSolutionSummary);}
                                        cisVal1counter++;
                                        }
                                    else{//make & use threads
                                        std::vector<std::thread> threadList;
                                        std::vector<int> solutionSameOrBetterTF;
                                        std::vector<FitnessMaximumSolutionSet> newSolutionSummaries;
                                        std::vector<SimplestRegPathIndividual> focalIndivsToTest;
                                        std::vector<genotypeSettings> gtypeSettingsList;
                                        std::vector<simulationSettings> simSettingsList;
                                        for(int t=0;t<threadsToUse;++t){//initialize
                                            solutionSameOrBetterTF.push_back(false);
                                            FitnessMaximumSolutionSet fmss;
                                            newSolutionSummaries.push_back(fmss);
                                            gtypeSet.cisVal1_=cisVal1counter;
                                            gtypeSettingsList.push_back(gtypeSet);
                                            simSettingsList.push_back(simSet);
                                            focalIndiv.SetGenotype(2,1,cisVal1counter);
                                            focalIndivsToTest.push_back(focalIndiv);
                                            cisVal1counter++;
                                            }//initialize
                                        for(int t=0;t<threadsToUse;++t){//make and run threads
                                            SimplestRegPathIndividual *fIndiv = &(focalIndivsToTest[t]);
                                            FitnessMaximumSolutionSet *fmss = &(newSolutionSummaries[t]);
                                            genotypeSettings *gset = &(gtypeSettingsList[t]);
                                            simulationSettings *simset = &(simSettingsList[t]);
                                            int *sameOrBetter = &(solutionSameOrBetterTF[t]);
                                            std::thread th(
                                                MaximizeUsingBitstringsOneReferenceGtype,fIndiv,simset,gset,fmss,
                                                                maxPopMeanFitness,sameOrBetter);
                                            threadList.push_back(move(th));
                                            }//t
                                        for(int t=0;t<threadList.size();++t){//synchronize threads
                                            threadList[t].join();}
                                        for(int t=0;t<threadList.size();++t){//collect data
                                            if(solutionSameOrBetterTF[t]){
                                                summariesOfSolutions.AddSolution(newSolutionSummaries[t]);
                                                maxPopMeanFitness=MAX(maxPopMeanFitness,newSolutionSummaries[t].wBarMax());
                                                }//solutionSameOrBetterTF
                                            }//t
                                        }//make & use threads
                                }//while cis1ValCounter
                            
                        }//cisVal0
                    //could end a thread here that maximized fitness over the whole cis parameter space for this TF pair
                    //each thread would end with a local FitnessMaximaSolutionSets object that stored all 2-allele cis genotypes having that maximum
                    //these would be concatenated to the global FitnessMaximaSolutionSets
                    }//tfVal1
                }//dosageVal1
            }//tfVal0
        }//dosageVal0
    coutLock.lock(); std::cout<<std::endl<<"**************************"<<std::endl; coutLock.unlock();
    summariesOfSolutions.PrintDataByPopt(outputfileSolutionSummaries, Popt);
    }//maximizeUsingBitstringsThreadable


void MaximizeUsingBitstringsThreadableAllCis(int bitstringLen, long double NtfsatPerAllele, long double deltaG1dosage,
                long double deltaG1, long double minExpression, long double maxExpression,
                long double Popt, long double omega,
                FitnessMaximaBitstringSolutions& wBarMaxAllSolutions, std::ostream& outputfileAllSolutions, bool saveAllSolutions,
                FitnessMaximaSolutionSets& summariesOfSolutions, std::ostream& outputfileSolutionSummaries,
                bool splitSinglePoptRun, uint64_t lowTf0dosage, uint64_t highTf0dosage){
    //this version puts the whole cis array into one thread
    coutLock.lock();
    std::cout<<"maximizing for Popt="<<Popt<<" and omega="<<omega;
    if(splitSinglePoptRun){
        std::cout<<" for tf0 dosage "<<lowTf0dosage;
        if(highTf0dosage>lowTf0dosage){
            std::cout<<" to "<<highTf0dosage;}
        }
    std::cout<<std::endl; coutLock.unlock();
    uint64_t maxBitstringVal = uint64_t(pow(2,bitstringLen));
     if(!splitSinglePoptRun){//to eliminate any ambiguity
        lowTf0dosage=0; highTf0dosage=maxBitstringVal-1;}
   
    simulationSettings simSet(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression,Popt,omega,
            splitSinglePoptRun,lowTf0dosage,highTf0dosage);
    genotypeSettings gtypeSet;
    SimplestRegPathIndividual focalIndiv(false);//indivAaBb
    long double maxPopMeanFitness=-one;
//    std::vector<long double> phat(10,-one), qhat(10,-one),popMeanPhenotypes(10,-one);
//    std::vector<bool> pNeutral(4,false), qNeutral(4,false);
    
//    for(uint64_t dosageVal0=0;dosageVal0<maxBitstringVal;++dosageVal0){//1st TF allele dosage  //dosageVal0=3
    for(uint64_t dosageVal0=lowTf0dosage;dosageVal0<=highTf0dosage;++dosageVal0){//1st TF allele dosage
        gtypeSet.dosageVal0_=dosageVal0;
        focalIndiv.SetGenotype(0,0,dosageVal0);
        for(uint64_t tfVal0=0;tfVal0<maxBitstringVal;++tfVal0){//1st TF allele product
            gtypeSet.tfVal0_=tfVal0;
            focalIndiv.SetGenotype(1,0,tfVal0);
            uint64_t highDosageVal1=dosageVal0;
            if(splitSinglePoptRun){//if so, will need to check every Tf1 dosage
                highDosageVal1=maxBitstringVal-1;}
//            for(uint64_t dosageVal1=0;dosageVal1<=dosageVal0;++dosageVal1){//2nd TF allele dosage // dosageVal1=3
            for(uint64_t dosageVal1=0;dosageVal1<=highDosageVal1;++dosageVal1){//2nd TF allele dosage
                gtypeSet.dosageVal1_=dosageVal1;
                focalIndiv.SetGenotype(0,1,dosageVal1);
                std::thread refThread;
                uint64_t tfVal1counter=0;
                while(tfVal1counter<=tfVal0){//tfVal1Counter=10
                    int availableThreads = MAX((unsigned int)1,refThread.hardware_concurrency());
                    int threadsToUse = MIN(availableThreads,(int)(tfVal0-tfVal1counter)+1);//+1 for when tfVal0==tfVal1
                    threadsToUse = MIN(threadsToUse,16);
                    if(threadsToUse<=1){//just call the function without using a thread
                            gtypeSet.tfVal1_=tfVal1counter;
                            focalIndiv.SetGenotype(1,1,tfVal1counter);
                            int solutionEqualsMaxOrBetter=false;
                            FitnessMaximaSolutionSets newSolutionSummaries(splitSinglePoptRun,lowTf0dosage,highTf0dosage);
                            MaximizeUsingBitstringsAllCisGtypes(&focalIndiv,&simSet,&gtypeSet,&newSolutionSummaries,
                                    maxPopMeanFitness,&solutionEqualsMaxOrBetter);
                            if(solutionEqualsMaxOrBetter){
                                maxPopMeanFitness=newSolutionSummaries.wBarMax(Popt);
                                summariesOfSolutions.ConcatenateSolutions(newSolutionSummaries);}
                            tfVal1counter++;
                            }
                        else{//make & use threads
                            std::vector<std::thread> threadList;
                            std::vector<int> solutionSameOrBetterTF;
                            std::vector<FitnessMaximaSolutionSets> newSolutionSummaries;
                            std::vector<SimplestRegPathIndividual> focalIndivsToTest;
                            std::vector<genotypeSettings> gtypeSettingsList;
                            std::vector<simulationSettings> simSettingsList;
                            for(int t=0;t<threadsToUse;++t){//initialize
                                solutionSameOrBetterTF.push_back(false);
                                FitnessMaximaSolutionSets fmss(splitSinglePoptRun,lowTf0dosage,highTf0dosage);
                                newSolutionSummaries.push_back(fmss);
                                gtypeSet.tfVal1_=tfVal1counter;
                                gtypeSettingsList.push_back(gtypeSet);
                                simSettingsList.push_back(simSet);
                                focalIndiv.SetGenotype(1,1,tfVal1counter);
                                focalIndivsToTest.push_back(focalIndiv);
                                tfVal1counter++;
                                }//initialize
                            for(int t=0;t<threadsToUse;++t){//make and run threads
                                SimplestRegPathIndividual *fIndiv = &(focalIndivsToTest[t]);
                                FitnessMaximaSolutionSets *fmss = &(newSolutionSummaries[t]);
                                genotypeSettings *gset = &(gtypeSettingsList[t]);
                                simulationSettings *simset = &(simSettingsList[t]);
                                int *sameOrBetter = &(solutionSameOrBetterTF[t]);
                                std::thread th(
                                    MaximizeUsingBitstringsAllCisGtypes,fIndiv,simset,gset,fmss,
                                                    maxPopMeanFitness,sameOrBetter);
                                threadList.push_back(move(th));
                                }//t
                            for(int t=0;t<threadList.size();++t){//synchronize threads
                                threadList[t].join();}//t=6
                            for(int t=0;t<threadList.size();++t){//collect data
                                if(solutionSameOrBetterTF[t]){
                                    if(summariesOfSolutions.wBarMax(Popt)<=newSolutionSummaries[t].wBarMax(Popt)){
                                        //probably not all 'same or better' threads will be equal, so skip those
                                        summariesOfSolutions.ConcatenateSolutions(newSolutionSummaries[t]);
                                        maxPopMeanFitness=MAX(maxPopMeanFitness,newSolutionSummaries[t].wBarMax(Popt));
                                        }
                                    }//solutionSameOrBetterTF
                                }//t
                            }//make & use threads
                    }//while tfVal0Counter
              }//dosageVal1
            }//tfVal0
        }//dosageVal0
    coutLock.lock(); std::cout<<std::endl<<"**************************"<<std::endl; coutLock.unlock();
    summariesOfSolutions.PrintDataByPopt(outputfileSolutionSummaries, Popt);
    }//MaximizeUsingBitstringsThreadableAllCis




void MaximizeUsingBitstrings(int bitstringLen, long double NtfsatPerAllele, long double deltaG1dosage,
                long double deltaG1, long double minExpression, long double maxExpression,
                long double Popt, long double omega,
                FitnessMaximaBitstringSolutions& wBarMaxAllSolutions, std::ostream& outputfileAllSolutions, bool saveAllSolutions,
                FitnessMaximaSolutionSets& summariesOfSolutions, std::ostream& outputfileSolutionSummaries,
                bool splitSinglePoptRun, uint64_t lowTf0dosage, uint64_t highTf0dosage){
    coutLock.lock();
    std::cout<<"maximizing for Popt="<<Popt<<" and omega="<<omega;
    if(splitSinglePoptRun){
        std::cout<<" for tf0 dosage "<<lowTf0dosage;
        if(highTf0dosage>lowTf0dosage){
            std::cout<<" to "<<highTf0dosage;}
        }
    std::cout<<std::endl; coutLock.unlock();
    uint64_t maxBitstringVal = uint64_t(pow(2,bitstringLen));
     if(!splitSinglePoptRun){//to eliminate any ambiguity
        lowTf0dosage=0; highTf0dosage=maxBitstringVal-1;}
    SimplestRegPathIndividual indivAABB(false), indivAABb(false), indivAAbb(false);//recombinants
    SimplestRegPathIndividual indivAaBB(false), indivAaBb(false), indivAabb(false);//
    SimplestRegPathIndividual indivaaBB(false), indivaaBb(false), indivaabb(false);
    SimplestRegPathIndividual focalIndiv(false);//indivAaBb
    long double wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,wBarMax;
    long double phAABB,phAABb,phAAbb,phAaBB,phAaBb,phAabb,phaaBB,phaaBb,phaabb;
    long double popMeanFitness,maxPopMeanFitness=-one;
    wAABB=wAABb=wAAbb=wAaBB=wAaBb=wAabb=waaBB=waaBb=waabb=wBarMax=popMeanFitness=zero;//fitnesses
    phAABB=phAABb=phAAbb=phAaBB=phAaBb=phAabb=phaaBB=phaaBb=phaabb=zero;//phenotypes
    std::vector<long double> phat(10,-one), qhat(10,-one),popMeanPhenotypes(10,-one);
    std::vector<bool> pNeutral(4,false), qNeutral(4,false);
    bool newBest=false;
    
    for(uint64_t dosageVal0=lowTf0dosage;dosageVal0<=highTf0dosage;++dosageVal0){//1st TF allele dosage
        focalIndiv.SetGenotype(0,0,dosageVal0);
        for(uint64_t tfVal0=0;tfVal0<maxBitstringVal;++tfVal0){//1st TF allele product
            focalIndiv.SetGenotype(1,0,tfVal0);
            uint64_t highDosageVal1=dosageVal0;
            if(splitSinglePoptRun){//if so, still need to check every Tf1 dosage
                highDosageVal1=maxBitstringVal-1;}
            for(uint64_t dosageVal1=0;dosageVal1<=highDosageVal1;++dosageVal1){//2nd TF allele dosage
                focalIndiv.SetGenotype(0,1,dosageVal1);
                for(uint64_t tfVal1=0;tfVal1<=tfVal0;++tfVal1){//2nd TF allele product
                    focalIndiv.SetGenotype(1,1,tfVal1);
                    //could put here a thread that maximized fitness over the whole cis parameter space for this TF pair
                    //so, if using threads, create a local FitnessMaximaSolutionSets object
                    
                    for(uint64_t cisVal0=0;cisVal0<maxBitstringVal;++cisVal0){//1st cis allele promoter
                        focalIndiv.SetGenotype(2,0,cisVal0);
                        for(uint64_t cisVal1=0;cisVal1<=cisVal0;++cisVal1){//2nd cis allele promoter
                            focalIndiv.SetGenotype(2,1,cisVal1);
                            for(int i=0;i<4;++i){pNeutral[i]=qNeutral[i]=false;}
                            wAABB=wAABb=wAAbb=wAaBB=wAaBb=wAabb=waaBB=waaBb=waabb=-one;
                            phAABB=phAABb=phAAbb=phAaBB=phAaBb=phAabb=phaaBB=phaaBb=phaabb=-one;
                            focalIndiv.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);
                            focalIndiv.CalculateFitness(Popt,omega);
                            std::string focalGtypeBitstringStr=focalIndiv.gtypeString(bitstringLen);
                            std::string focalGtypeMismatchStr=focalIndiv.mismatchStringMathematicaFormat();
                            //create recombinant genotypes
                            int numMaxima=0;
                            if(focalIndiv.IsTFheterozygote() && focalIndiv.IsCisHeterozygote()){//maximize for p & q
                                    indivAaBb=focalIndiv;
                                    
                                    indivAABB.TFdosage_[0]=indivAABB.TFdosage_[1]=dosageVal0;
                                    indivAABB.TFproduct_[0]=indivAABB.TFproduct_[1]=tfVal0;
                                    indivAABB.cis_[0]=indivAABB.cis_[1]=cisVal0;
                                
                                    indivAABb.TFdosage_[0]=indivAABb.TFdosage_[1]=dosageVal0;
                                    indivAABb.TFproduct_[0]=indivAABb.TFproduct_[1]=tfVal0;
                                    indivAABb.cis_[0]=cisVal0; indivAABb.cis_[1]=cisVal1;
                                
                                    indivAAbb.TFdosage_[0]=indivAAbb.TFdosage_[1]=dosageVal0;
                                    indivAAbb.TFproduct_[0]=indivAAbb.TFproduct_[1]=tfVal0;
                                    indivAAbb.cis_[0]=indivAAbb.cis_[1]=cisVal1;
                                
                                    indivAaBB.TFdosage_[0]=dosageVal0; indivAaBB.TFdosage_[1]=dosageVal1;
                                    indivAaBB.TFproduct_[0]=tfVal0; indivAaBB.TFproduct_[1]=tfVal1;
                                    indivAaBB.cis_[0]=indivAaBB.cis_[1]=cisVal0;
                                
                                    indivAabb.TFdosage_[0]=dosageVal0; indivAabb.TFdosage_[1]=dosageVal1;
                                    indivAabb.TFproduct_[0]=tfVal0; indivAabb.TFproduct_[1]=tfVal1;
                                    indivAabb.cis_[0]=indivAabb.cis_[1]=cisVal1;
                                
                                    indivaaBB.TFdosage_[0]=indivaaBB.TFdosage_[1]=dosageVal1;
                                    indivaaBB.TFproduct_[0]=indivaaBB.TFproduct_[1]=tfVal1;
                                    indivaaBB.cis_[0]=indivaaBB.cis_[1]=cisVal0;
                                
                                    indivaaBb.TFdosage_[0]=indivaaBb.TFdosage_[1]=dosageVal1;
                                    indivaaBb.TFproduct_[0]=indivaaBb.TFproduct_[1]=tfVal1;
                                    indivaaBb.cis_[0]=cisVal0; indivaaBb.cis_[1]=cisVal1;
                                
                                    indivaabb.TFdosage_[0]=indivaabb.TFdosage_[1]=dosageVal1;
                                    indivaabb.TFproduct_[0]=indivaabb.TFproduct_[1]=tfVal1;
                                    indivaabb.cis_[0]=indivaabb.cis_[1]=cisVal1;
    

                                    indivAABB.Reset();indivAABb.Reset();indivAAbb.Reset();indivAaBB.Reset();
                                    indivAabb.Reset();indivaaBB.Reset();indivaaBb.Reset();indivaabb.Reset();
                                    std::string focal=focalIndiv.mismatchStringMathematicaFormat();
                                    std::string mAABB=indivAABB.mismatchStringMathematicaFormat();
                                    std::string mAABb=indivAABb.mismatchStringMathematicaFormat();
                                    std::string mAAbb=indivAAbb.mismatchStringMathematicaFormat();
                                    std::string mAaBB=indivAaBB.mismatchStringMathematicaFormat();
                                    std::string mAaBb=focalIndiv.mismatchStringMathematicaFormat();
                                    std::string mAabb=indivAabb.mismatchStringMathematicaFormat();
                                    std::string maaBB=indivaaBB.mismatchStringMathematicaFormat();
                                    std::string maaBb=indivaaBb.mismatchStringMathematicaFormat();
                                    std::string maabb=indivaabb.mismatchStringMathematicaFormat();
                                    indivAABB.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);
                                    indivAABb.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);
                                    indivAAbb.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);
                                    indivAaBB.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);
                                    indivAabb.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);
                                    indivaaBB.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);
                                    indivaaBb.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);
                                    indivaabb.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);

                                    phAABB=indivAABB.phenotype();
                                    phAABb=indivAABb.phenotype();
                                    phAAbb=indivAAbb.phenotype();
                                    phAaBB=indivAaBB.phenotype();
                                    phAaBb=focalIndiv.phenotype();
                                    phAabb=indivAabb.phenotype();
                                    phaaBB=indivaaBB.phenotype();
                                    phaaBb=indivaaBb.phenotype();
                                    phaabb=indivaabb.phenotype();

                                    wAABB=indivAABB.CalculateFitness(Popt,omega);
                                    wAABb=indivAABb.CalculateFitness(Popt,omega);
                                    wAAbb=indivAAbb.CalculateFitness(Popt,omega);
                                    wAaBB=indivAaBB.CalculateFitness(Popt,omega);
                                    wAaBb=focalIndiv.fitness();
                                    wAabb=indivAabb.CalculateFitness(Popt,omega);
                                    waaBB=indivaaBB.CalculateFitness(Popt,omega);
                                    waaBb=indivaaBb.CalculateFitness(Popt,omega);
                                    waabb=indivaabb.CalculateFitness(Popt,omega);
                                
                                //maximize for p & q
                                MaximizePopMeanFitnessPandQv2(wAABB,wAABb,wAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,popMeanFitness,phat,qhat,pNeutral,qNeutral,numMaxima,focalGtypeBitstringStr);
                                collectMeanPhenotypesPQ(phAABB,phAABb,phAAbb,phAaBB,phAaBb,phAabb,phaaBB,phaaBb,phaabb,phat,qhat,popMeanPhenotypes,numMaxima);
                                    }
                                else if(focalIndiv.IsTFheterozygote()){//maximize for p
                                    indivAaBB=focalIndiv;
                                    
                                    indivAABB.TFdosage_[0]=indivAABB.TFdosage_[1]=dosageVal0;
                                    indivAABB.TFproduct_[0]=indivAABB.TFproduct_[1]=tfVal0;
                                    indivAABB.cis_[0]=indivAABB.cis_[1]=cisVal0;

                                    indivaaBB.TFdosage_[0]=indivaaBB.TFdosage_[1]=dosageVal1;
                                    indivaaBB.TFproduct_[0]=indivaaBB.TFproduct_[1]=tfVal1;
                                    indivaaBB.cis_[0]=indivaaBB.cis_[1]=cisVal0;
                                    
                                    indivAABB.Reset();indivaaBB.Reset();
                                    indivAABB.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);
                                    indivaaBB.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);

                                    phAABB=indivAABB.phenotype();
                                    phAaBB=focalIndiv.phenotype();
                                    phaaBB=indivaaBB.phenotype();

                                    wAABB=indivAABB.CalculateFitness(Popt,omega);
                                    wAaBB=focalIndiv.fitness();
                                    waaBB=indivaaBB.CalculateFitness(Popt,omega);
                                    MaximizePopMeanFitnessP(wAABB,wAaBB,waaBB,popMeanFitness,phat,pNeutral,numMaxima);
                                    for(int m=0;m<numMaxima;++m){qhat[m]=one;qNeutral[m]=false;}
                                    collectMeanPhenotypesP(phAABB,phAaBB,phaaBB,phat,popMeanPhenotypes,numMaxima);
                                    }
                                else if(focalIndiv.IsCisHeterozygote()){//maximize for q
                                    indivAABb = focalIndiv;
                                    
                                    indivAABB.TFdosage_[0]=indivAABB.TFdosage_[1]=dosageVal0;
                                    indivAABB.TFproduct_[0]=indivAABB.TFproduct_[1]=tfVal0;
                                    indivAABB.cis_[0]=indivAABB.cis_[1]=cisVal0;

                                    indivAAbb.TFdosage_[0]=indivAAbb.TFdosage_[1]=dosageVal0;
                                    indivAAbb.TFproduct_[0]=indivAAbb.TFproduct_[1]=tfVal0;
                                    indivAAbb.cis_[0]=indivAAbb.cis_[1]=cisVal1;

                                    indivAABB.Reset();indivAAbb.Reset();
                                    indivAABB.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);
                                    indivAAbb.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);
                                    phAABB=indivAABB.phenotype();
                                    phAABb=focalIndiv.phenotype();
                                    phAAbb=indivAAbb.phenotype();

                                    wAABB=indivAABB.CalculateFitness(Popt,omega);
                                    wAABb=focalIndiv.fitness();
                                    wAAbb=indivAAbb.CalculateFitness(Popt,omega);
                                    MaximizePopMeanFitnessQ(wAABB,wAABb,wAAbb,popMeanFitness,qhat,qNeutral,numMaxima);
                                    for(int m=0;m<numMaxima;++m){phat[m]=one;pNeutral[m]=false;}
                                    collectMeanPhenotypesQ(phAABB,phAABb,phAAbb,qhat,popMeanPhenotypes,numMaxima);
                                    }
                                else{//double homozygote AABB
                                    indivAABB=focalIndiv;
                                    phAABB=focalIndiv.phenotype();
                                    wAABB=focalIndiv.fitness();
                                    phat[0]=qhat[0]=one;
                                    pNeutral[0]=qNeutral[0]=false;
                                    numMaxima=1;
                                    popMeanFitness=focalIndiv.fitness();
                                    popMeanPhenotypes[0]=focalIndiv.phenotype();
                                    }

                            popMeanFitness=MIN(MAX(ROUND(popMeanFitness,decimalDigitsToRound),zero),one);
                            for(int m=0;m<numMaxima;++m){
                                popMeanPhenotypes[m]=MIN(MAX(ROUND(popMeanPhenotypes[m],decimalDigitsToRound),zero),one);}

                            
                            
                            
                            if(popMeanFitness==maxPopMeanFitness){

                                for(int m=0;m<numMaxima;++m){
                                    if(phat[m]==-one && qhat[m]==-one) break;//reached the end
                                    std::string ht("___"), phatstr, qhatstr;//, ABgtype;
                                    std::stringstream phatstrSS,qhatstrSS;
                                    phatstrSS<<phat[m]; qhatstrSS<<qhat[m];
//                                    phatstr=std::to_string(phat[m]); qhatstr=std::to_string(qhat[m]);
                                    phatstr=phatstrSS.str(); qhatstr=phatstrSS.str();
                                    SimplestRegPathIndividual solutionIndiv;
                                    BitstringGenotypeData referenceIndivData, solutionIndivData;
                                    if(phat[m]==one){//the A TF allele is fixed
                                            if(qhat[m]==one){//the B cis allele is fixed
                                                    solutionIndiv=indivAABB;
                                                    }
                                                else if(qhat[m]==zero){//the b cis allele is fixed
                                                    solutionIndiv=indivAAbb;
                                                    }
                                                else if(qNeutral[m]){
                                                    qhatstr="neutral";
                                                    ht[2]='n';
                                                    solutionIndiv=indivAABb;
                                                    }
                                                else{//cis locus is heterozygous
                                                    ht[2]='c';
                                                    solutionIndiv=indivAABb;
                                                    }
                                            }
                                        else if(phat[m]==zero){//the a TF allele is fixed
                                            if(qhat[m]==one){//the B cis allele is fixed
                                                    solutionIndiv=indivaaBB;
                                                    }
                                                else if(qhat[m]==zero){//the b cis allele is fixed
                                                    solutionIndiv=indivaabb;
                                                    }
                                                else if(qNeutral[m]){
                                                    qhatstr="neutral";
                                                    ht[2]='n';
                                                    solutionIndiv=indivaaBb;
                                                    }
                                                else{//cis locus is heterozygous
                                                    ht[2]='c';
                                                    solutionIndiv=indivaaBb;
                                                    }
                                            }
                                        else if(pNeutral[m]){
                                             phatstr="neutral";
                                            //figure out which piece of the TF is heterozygous, if either
                                            if(dosageVal0 != dosageVal1){
                                                ht[0]='n';}//dosages differ
                                            if(tfVal0 != tfVal1){
                                                ht[1]='n';}//product sites differ
                                            if(qhat[m]==one){//the B cis allele is fixed
                                                    solutionIndiv=indivAaBB;
                                                    }
                                                else if(qhat[m]==zero){//the b cis allele is fixed
                                                    solutionIndiv=indivAabb;
                                                    }
                                                else if(qNeutral[m]){
                                                    qhatstr="neutral";
                                                    ht[2]='n';
                                                    solutionIndiv=focalIndiv;//AaBb
                                                    }
                                                else{//cis locus is heterozygous
                                                    ht[2]='c';
                                                    solutionIndiv=focalIndiv;//AaBb
                                                    }
                                            }//neutral TF
                                        else{//the TF is polymorphic & non-neutral
                                            //figure out which pieces of the TF are heterozygous/neutral
                                            if(dosageVal0 != dosageVal1){
                                                ht[0]='d';}//dosages differ and aren't neutral
                                            if(qhat[m]==one){//the B cis allele is fixed
                                                    solutionIndiv=indivAaBB;
                                                    }
                                                else if(qhat[m]==zero){//the b cis allele is fixed
                                                    solutionIndiv=indivAabb;
                                                    }
                                                else if(qNeutral[m]){
                                                    ht[2]='n';
                                                    qhatstr="neutral";
                                                    solutionIndiv=focalIndiv;//AaBb
                                                    }
                                                else{//cis locus is heterozygous
                                                    ht[2]='c';
                                                    solutionIndiv=focalIndiv;//AaBb
                                                    }
                                            if(tfVal0 != tfVal1){//product sites differ
                                                solutionIndiv.SetMismatchesUsingBitstrings();
                                                if(solutionIndiv.IsTFprodHeterozygote(true)){//mismatch heterozygote
                                                        ht[1]='p';}
                                                    else{
                                                        ht[1]='n';}
                                                    }
                                            }//polymorphic TF

                                    referenceIndivData.CollectData(focalIndiv,bitstringLen);
                                    solutionIndiv.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);//sets mismatch pattern
                                    //see if this solution has been found via another reference sequence
                                    long double pMostCommon=-one,qMostCommon=-one;
                                    if(!pNeutral[m]){pMostCommon=MAX(phat[m],one-phat[m]);}
                                    if(!qNeutral[m]){qMostCommon=MAX(qhat[m],one-qhat[m]);}
                                    solutionIndiv.SetMismatchesUsingBitstrings();
                                    std::string solGtype=solutionIndiv.gtypeString(bitstringLen);
                                    std::string refGtype=focalIndiv.gtypeString(bitstringLen);
                                    std::string mhp = solutionIndiv.mismatchHetType();
                                    std::string mp = solutionIndiv.mismatchStringMathematicaFormat();
                                    int hc = codeToInt(ht), mhc = codeToInt(mhp);
                                    FitnessMaximumSolutionSet thisSolutionSummary(Popt,omega,NtfsatPerAllele,bitstringLen,popMeanFitness,
                                                                                  popMeanPhenotypes[m],pMostCommon,
                                                                                  qMostCommon,pNeutral[m],qNeutral[m],
                                                                                  hc,ht,mhc,mhp,mp,solGtype,refGtype,
                                                                                  splitSinglePoptRun,lowTf0dosage,highTf0dosage);
                                    summariesOfSolutions.AddSolution(thisSolutionSummary);
                                    }//m
                                }//popMeanFitness==maxPopMeanFitness
 
                            if(popMeanFitness>maxPopMeanFitness){
                                newBest=true;
                                for(int m=0;m<numMaxima;++m){
                                    SimplestRegPathIndividual solutionIndiv;
                                    BitstringGenotypeData referenceIndivData, solutionIndivData;

                                    std::string ht("___"), phatstr, qhatstr;
                                    phatstr=std::to_string(phat[m]); qhatstr=std::to_string(qhat[m]);
                                    if(phat[m]==one){//the A TF allele is fixed
                                            if(qhat[m]==one){//the B cis allele is fixed
                                                    solutionIndiv=indivAABB;
                                                    }
                                                else if(qhat[m]==zero){//the b cis allele is fixed
                                                    solutionIndiv=indivAAbb;
                                                    }
                                                else if(qNeutral[m]){
                                                    qhatstr="neutral";
                                                    ht[2]='n';
                                                    solutionIndiv=indivAABb;
                                                    }
                                                else{//cis locus is heterozygous
                                                    ht[2]='c';
                                                    solutionIndiv=indivAABb;
                                                    }
                                            }
                                        else if(phat[m]==zero){//the a TF allele is fixed
                                            if(qhat[m]==one){//the B cis allele is fixed
                                                    solutionIndiv=indivaaBB;
                                                    }
                                                else if(qhat[m]==zero){//the b cis allele is fixed
                                                    solutionIndiv=indivaabb;
                                                    }
                                                else if(qNeutral[m]){
                                                    qhatstr="neutral";
                                                    ht[2]='n';
                                                    solutionIndiv=indivaaBb;
                                                    }
                                                else{//cis locus is heterozygous
                                                    ht[2]='c';
                                                    solutionIndiv=indivaaBb;
                                                    }
                                            }
                                        else if(pNeutral[m]){
                                            phatstr="neutral";
                                            //figure out which piece of the TF is heterozygous, if either
                                            if(dosageVal0 != dosageVal1){
                                                ht[0]='n';}//dosages differ; neutrality is here
                                            if(tfVal0 != tfVal1){
                                                ht[1]='n';}//product sites differ; neutrality is here
                                            if(qhat[m]==one){//the B cis allele is fixed
                                                    solutionIndiv=indivAaBB;
                                                    }
                                                else if(qhat[m]==zero){//the b cis allele is fixed
                                                    solutionIndiv=indivAabb;
                                                    }
                                                else if(qNeutral[m]){//shouldn't get here
                                                    qhatstr="neutral";
                                                    ht[2]='n';
                                                    solutionIndiv=focalIndiv;
                                                    }
                                                else{//cis locus is heterozygous
                                                    ht[2]='c';
                                                    solutionIndiv=focalIndiv;
                                                    }
                                            }//neutral TF
                                        else{//the TF is polymorphic & non-neutral
                                            //figure out which pieces of the TF are heterozygous/neutral
                                            if(dosageVal0 != dosageVal1){
                                                ht[0]='d';}//dosages differ and aren't neutral
                                            if(qhat[m]==one){//the B cis allele is fixed
                                                    solutionIndiv=indivAaBB;
                                                    }
                                                else if(qhat[m]==zero){//the b cis allele is fixed
                                                    solutionIndiv=indivAabb;
                                                    }
                                                else if(qNeutral[m]){
                                                    ht[2]='n';
                                                    qhatstr="neutral";
                                                    solutionIndiv=focalIndiv;//AaBb
                                                    }
                                                else{//cis locus is heterozygous
                                                    ht[2]='c';
                                                    solutionIndiv=focalIndiv;//AaBb
                                                    }
                                            if(tfVal0 != tfVal1){//product sites differ
                                                solutionIndiv.SetMismatchesUsingBitstrings();
                                                if(solutionIndiv.IsTFprodHeterozygote(true)){//mismatch heterozygote
                                                        ht[1]='p';}
                                                    else{
                                                        ht[1]='n';}
                                                    }
                                            }//polymorphic TF
                                    maxPopMeanFitness=popMeanFitness;
                                    referenceIndivData.CollectData(focalIndiv,bitstringLen);
                                    solutionIndiv.CalculatePhenotype(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression);//sets mismatch pattern data
                                    //see if this solution has been found via another reference sequence
                                    long double pMostCommon=-one,qMostCommon=-one;
                                    if(!pNeutral[m]){pMostCommon=MAX(phat[m],one-phat[m]);}
                                    if(!qNeutral[m]){qMostCommon=MAX(qhat[m],one-qhat[m]);}
                                    solutionIndiv.SetMismatchesUsingBitstrings();
                                    std::string solGtype=solutionIndiv.gtypeString(bitstringLen);
                                    std::string refGtype=focalIndiv.gtypeString(bitstringLen);
                                    std::string mhp = solutionIndiv.mismatchHetType();
                                    std::string mp = solutionIndiv.mismatchStringMathematicaFormat();
                                    int hc = codeToInt(ht), mhc = codeToInt(mhp);
                                    FitnessMaximumSolutionSet thisSolutionSummary(Popt,omega,NtfsatPerAllele,bitstringLen,popMeanFitness,
                                                                                  popMeanPhenotypes[m],pMostCommon,
                                                                                  qMostCommon,pNeutral[m],qNeutral[m],
                                                                                  hc,ht,mhc,mhp,mp,solGtype,refGtype,
                                                                                  splitSinglePoptRun,lowTf0dosage,highTf0dosage);
                                    summariesOfSolutions.AddSolution(thisSolutionSummary);
                                    }//m
                                }//popMeanFitness>maxPopMeanFitness
                            }//cisVal1
                        }//cisVal0
                    //could end a thread here that maximized fitness over the whole cis parameter space for this TF pair
                    //each thread would end with a local FitnessMaximaSolutionSets object that stored all 2-allele cis genotypes having that maximum
                    //these would be concatenated to the global FitnessMaximaSolutionSets
                    }//tfVal1
                }//dosageVal1
            }//tfVal0
        }//dosageVal0
    coutLock.lock(); std::cout<<std::endl<<"**************************"<<std::endl; coutLock.unlock();
    summariesOfSolutions.PrintDataByPopt(outputfileSolutionSummaries, Popt);
    }//maximizeUsingBitstrings






std::string binaryFromUnsignedLongLong(unsigned long long val, bool fullString, int prependZerosToLen){
    int i,last1=0;
    size_t length=CHAR_BIT*sizeof(unsigned long long);
    std::string bitstring(length,'0');
    unsigned int lastbitMask=1;
    for(i=0;i<length;++i){
        bitstring[length-1-i] ='0'+char(val & lastbitMask);
        if(bitstring[length-1-i]=='1') {last1=i+1;}
        val = val>>1;
        }
    if(fullString) return bitstring;
    bitstring.erase(bitstring.begin(),bitstring.end()-last1);
    std::string zeros(prependZerosToLen-last1,'0');
    zeros+=bitstring;
    return zeros;
    }//binaryFromUnsignedLongLong


void groupGtypesByMismatches(std::vector<std::string>& allGtypes, std::vector<std::string>& allMismatches,
                             std::vector<std::string>& allHetStatuses,
                             std::vector<std::string>& uniqueMismatches, std::vector<std::vector<std::string> >& gtypesPerMismatch,
                             std::vector<std::vector<std::string> >& hetStatusPerMismatch){
    uniqueMismatches.clear();
    gtypesPerMismatch.clear();
    hetStatusPerMismatch.clear();
    std::string gtype,mismatch,hetstatus;
    bool found;
    for(unsigned long m=0;m<allMismatches.size();++m){
        found=false;
        gtype=allGtypes[m];
        mismatch=allMismatches[m];
        hetstatus=allHetStatuses[m];
        for(unsigned long u=0;u<uniqueMismatches.size();++u){
            std::string testMismatch(uniqueMismatches[u]);
            if(mismatch==testMismatch){
                gtypesPerMismatch[u].push_back(gtype);
                hetStatusPerMismatch[u].push_back(hetstatus);
                found=true;
                break;}
            }//u
        if(found){continue;}
        //make a new place for it
        uniqueMismatches.push_back(mismatch);
        std::vector<std::string> gtypeSet, hetStatusSet;
        gtypeSet.push_back(gtype);
        gtypesPerMismatch.push_back(gtypeSet);
        hetStatusSet.push_back(hetstatus);
        hetStatusPerMismatch.push_back(hetStatusSet);
        }//m
    }//groupGtypesByMismatches

/*
void DetermineViableMismatchCombinations(int bitstringLen){
    //this cycles through bitstrings and calculates the mismatch universe
    std::string tf1,tf2,cis1,cis2, m11,m21,m12,m22,tfhet,cishet,hetstatus,gtype,mismatches;
    int mTf1cis1,mTf1cis2,mTf2cis1,mTf2cis2;//mismatches
    bool tfHomozygote, cisHomozygote;
    std::vector<std::string> allGtypes, allMismatches, uniqueMismatches, allHetStatuses;
    std::vector<std::vector<std::string> > gtypesPerMismatch,hetStatusPerMismatch;
    uint64_t maxBitVal = uint64_t(pow(2,bitstringLen));
    for(uint64_t tf1val=0;tf1val<maxBitVal;++tf1val){
        tf1=binaryFromUnsignedLongLong(tf1val,false,3);
        for(uint64_t tf2val=tf1val;tf2val<maxBitVal;++tf2val){
            tf2=binaryFromUnsignedLongLong(tf2val,false,3);
            tfHomozygote = (tf1val==tf2val);
            if(tfHomozygote){tfhet="_";} else{tfhet="t";}
            for(uint64_t cis1val=0;cis1val<maxBitVal;++cis1val){
                cis1=binaryFromUnsignedLongLong(cis1val,false,3);
                mTf1cis1= (int)HammingDistance(tf1val,cis1val);
                mTf2cis1= (int)HammingDistance(tf2val,cis1val);
                m11=std::to_string(mTf1cis1);
                m21=std::to_string(mTf2cis1);
                for(uint64_t cis2val=cis1val;cis2val<maxBitVal;++cis2val){
                    cis2=binaryFromUnsignedLongLong(cis2val,false,3);
                    cisHomozygote = (cis1val==cis2val);
                    if(cisHomozygote){cishet="_";} else{cishet="c";}
                    hetstatus=tfhet+cishet;
                    mTf1cis2= (int)HammingDistance(tf1val,cis2val);
                    mTf2cis2= (int)HammingDistance(tf2val,cis2val);
                    m12=std::to_string(mTf1cis2);
                    m22=std::to_string(mTf2cis2);
                    gtype="[["+tf1+","+tf2+"],["+cis1+","+cis2+"]]";
                    mismatches="{{"+m11+","+m21+"},{"+m12+","+m22+"}}";
                    allGtypes.push_back(gtype);
                    allMismatches.push_back(mismatches);
                    allHetStatuses.push_back(hetstatus);
                    }//cis2val
                }//cis1val
            }//tf2val
        }//tf1val
    groupGtypesByMismatches(allGtypes,allMismatches,allHetStatuses,uniqueMismatches,gtypesPerMismatch,hetStatusPerMismatch);
    for(unsigned long m=0;m<uniqueMismatches.size();++m){
        std::cout<<uniqueMismatches[m]<<std::endl;
        for(unsigned long g=0;g<gtypesPerMismatch[m].size();++g){
            std::cout<<'\t'<<gtypesPerMismatch[m][g]<<'\t'<<hetStatusPerMismatch[m][g]<<std::endl;}
        }//m
    std::cout<<std::endl;
    }//DetermineViableMismatchCombinations
*/


void DetermineViableMismatchCombinations(int bitstringLen){
    //this cycles through bitstrings and calculates the mismatch universe, using stringstreams instead of strings
    std::string tf1,tf2,cis1,cis2,tfhet,cishet,hetstatus,gtype,mismatches;//m11,m21,m12,m22
    int mTf1cis1,mTf1cis2,mTf2cis1,mTf2cis2;//mismatches
    bool tfHomozygote, cisHomozygote;
    std::vector<std::string> allGtypes, allMismatches, uniqueMismatches, allHetStatuses;
    std::vector<std::vector<std::string> > gtypesPerMismatch,hetStatusPerMismatch;
    uint64_t maxBitVal = uint64_t(pow(2,bitstringLen));
    for(uint64_t tf1val=0;tf1val<maxBitVal;++tf1val){
        tf1=binaryFromUnsignedLongLong(tf1val,false,3);
        for(uint64_t tf2val=tf1val;tf2val<maxBitVal;++tf2val){
            tf2=binaryFromUnsignedLongLong(tf2val,false,3);
            tfHomozygote = (tf1val==tf2val);
            if(tfHomozygote){tfhet="_";} else{tfhet="t";}
            for(uint64_t cis1val=0;cis1val<maxBitVal;++cis1val){
                cis1=binaryFromUnsignedLongLong(cis1val,false,3);
                mTf1cis1= (int)HammingDistance(tf1val,cis1val);
                mTf2cis1= (int)HammingDistance(tf2val,cis1val);
                for(uint64_t cis2val=cis1val;cis2val<maxBitVal;++cis2val){
                    cis2=binaryFromUnsignedLongLong(cis2val,false,3);
                    cisHomozygote = (cis1val==cis2val);
                    if(cisHomozygote){cishet="_";} else{cishet="c";}
                    hetstatus=tfhet+cishet;
                    mTf1cis2= (int)HammingDistance(tf1val,cis2val);
                    mTf2cis2= (int)HammingDistance(tf2val,cis2val);
                    std::stringstream gtypeSS,mismatchesSS;
                    gtypeSS<<"[["<<tf1<<","<<tf2<<"],["<<cis1<<","<<cis2<<"]]";
                    mismatchesSS<<"{{"<<mTf1cis1<<","<<mTf2cis1<<"},{"<<mTf1cis2<<","<<mTf2cis2<<"}}";
                    gtype=gtypeSS.str();
                    mismatches=mismatchesSS.str();
                    allGtypes.push_back(gtype);
                    allMismatches.push_back(mismatches);
                    allHetStatuses.push_back(hetstatus);
                    }//cis2val
                }//cis1val
            }//tf2val
        }//tf1val
    groupGtypesByMismatches(allGtypes,allMismatches,allHetStatuses,uniqueMismatches,gtypesPerMismatch,hetStatusPerMismatch);
    coutLock.lock();
    for(unsigned long m=0;m<uniqueMismatches.size();++m){
        std::cout<<uniqueMismatches[m]<<std::endl;
        for(unsigned long g=0;g<gtypesPerMismatch[m].size();++g){
            std::cout<<'\t'<<gtypesPerMismatch[m][g]<<'\t'<<hetStatusPerMismatch[m][g]<<std::endl;}
        }//m
    std::cout<<std::endl;
    coutLock.unlock();
    }//DetermineViableMismatchCombinationsSS



using namespace std;
int main(int argc, const char * argv[]) {
    /*
    call in unix:
        $ g++ main.cpp -o fitnessOverdomOptGtype
        chmod fitnessOverdomOptGtype+x fitnessOverdomOptGtype.out
        ./fitnessOverdomOptGtype [bitstringLen] [PoptLow] [PoptHigh] [PoptSteps] [PoptStepSize]
        example:
        ./fitnessOverdomOptGtype 3 0 1000 1000 50
            #this runs the whole simulation with bitstringLen=3, starting at Popt=0, ending at Popt=1, in steps of 0.05
        ./fitnessOverdomOptGtype 3 10 20 1000 5
            #this runs the simulation with bitstringLen=3, but just the parts starting at Popt=0.01 and ending at Popt=0.02 in steps of 0.005
        ./fitnessOverdomOptGtype 3 10 10 1000 5
            #this runs the simulation with bitstringLen=3, but just for Popt=0.01 (steps and step size are irrelevant)
        ./fitnessOverdomOptGtype 3 1 1 100 1
            #this does the same thing
    */
        
/*
			//for debugging using cases from Mathematica file "wBar derivatives.nb"
			//plotFitnessSurface[wAABB_, wAABb_, wAAbb_, wAaBB_, wAaBb_, wAabb_,waaBB_, waaBb_, waabb_] :=
			//		Plot3D[wbar[wwAABB, wAABb, wAAbb, wAaBB, wAaBb, wAabb,waaBB, waaBb, waabb, p, q], {p, 0, 1}, {q, 0, 1}, PlotRange -> All];
			//wbar[wAABB_, wAABb_, wAAbb_, wAaBB_, wAaBb_, wAabb_,waaBB_, waaBb_, waabb_, p_, q_] :=
			//		wAABB p^2 q^2 + 2 wAABb p^2 q (1 - q) + wAAbb p^2 (1 - q)^2 +
			//		2 wAaBB  p (1 - p) q^2 + 4 wAaBb p (1 - p) q (1 - q) +
			//		2 wAabb p (1 - p) (1 - q)^2 + waaBB (1 - p)^2 q^2 +
			//		2 waaBb (1 - p)^2 q (1 - q) + waabb (1 - p)^2 (1 - q)^2;

	
			//plotFitnessSurface[0.5, 0.4, 0.6, 0.5, 0.2, 0.3, 0.7, 0.4, 0.5]//concave throughout; wbar->0.7,phat->1,qhat->0 --> correctly maximized
			//plotFitnessSurface[0.5, 0.6, 0.4, 0.5, 0.7, 0.65, 0.3, 0.35, 0.45]//convex throughout;
				//Mathematica values: wbar->0.558546,phat->0.324727,qhat->0.574640
				//close! maximized here to wbar->0.558546,phat->0.324729,qhat->0.574644
			//plotFitnessSurface2[0.5, 0.6, 0.4, 0.5, 0.5, 0.65, 0.3, 0.6, 0.45]//saddle shape
				//Mathematica finds the wrong local maximum: wbar->0.533333,p->0,q->0.333323
				//correctly maximized here to wbar->0.538889,phat->0.44444,qhat->0
			//plotFitnessSurface[0.05, 0.03, 0.06, 0.02, 0.9, 0.025, 0.06, 0.035, 0.065]
				//local maxima in the middle and in all 4 corners; middle highest
				//Mathematica values: wbar->0.253461, p -> 0.496251, q -> 0.49634
				//close! maximized here to wbar=0.253461, phat=0.496535, qhat=0.496378
			//plotFitnessSurface[0.25, 0.03, 0.06, 0.02, 1, 0.025, 0.06, 0.035, 0.065]
				//convex hump in the middle, but corner {1,1} is high at wbar=0.25
				//Mathematica values are wrong: wbar->0.25, phat -> 1., qhat -> 1.
				//correctly maximized here to wbar->0.6, {p->0, q->1} & {p->1, q->0}
				//plotFitnessSurface[0.5, 0.3, 0.6, 0.2, 1, 0.25, 0.6, 0.35, 0.65]
					//convex hump in the middle, but corner {0,0} is high at wbar=0.65
					//Mathematica values: wbar->0.65, phat -> 0, qhat -> 0
					//correctly maximized here to wbar->0.6, {p->0, q->1} & {p->1, q->0}
			//plotFitnessSurface[0.5, 0.3, 0.6, 0.2, 1, 0.25, 0.6, 0.35, 0.05]
				//convex hump in the middle, but corners {1,0} and {0,1} are high at wbar=0.6
				//Mathematica erroneously finds the central peak:
					//wbar->0.501564, phat -> 0.574541, qhat -> 0.534984
				//correctly maximized here to wbar->0.6, {p->0, q->1} & {p->1, q->0}

	//plotFitnessSurface[0.05, 0.6, 0.06, 0.02, 0.9, 0.025, 0.06, 0.035, 0.065]
		//convex hump in the middle, with one convex edge having slightly lower mean fitness
		//Mathematica correctly finds the central peak:
			//wbar->0.360261, phat -> 0.755171, qhat -> 0.496902
		//correctly maximized here to wbar->0.6, {p->0, q->1} & {p->1, q->0}
		
			//plotFitnessSurface[0.22, 0.13, 0.26, 0.2, 0.5, 0.24, 0.23, 0.15, 0.25]//convex hump in the middle
				//Mathematica values are wrong: wbar->0.25,p->1,q->1
				//correctly maximized here to wbar->0.276747, p->0.511748, q->0.596420

			//plotFitnessSurface2[0, 0, 0.0001389, 0, 0, 0.954332, 0, 0.0720475, 0]//wbar=0.477201, p=0.499964, q=1
   
            //plotFitnessSurface2[0.016, 0.546, 0.797, 0.958, 0.958, 0.958, 0.797, 0.546, 0.016]
                        //solutions: wBar=0.8205; {p,q}=>{0.145965,1} and (0.854035,0}
			long double wAABB=(long double)0.016;
			long double wAABb=(long double)0.546;
			long double WAAbb=(long double)0.797;
			long double wAaBB=(long double)0.958;
			long double wAaBb=(long double)0.958;
			long double wAabb=(long double)0.958;
			long double waaBB=(long double)0.797;
			long double waaBb=(long double)0.546;
			long double waabb=(long double)0.016;
			long double wBarMax=zero;
			std::vector<long double> phat,qhat;
            std::vector<bool> pNeutral,qNeutral;
			for(int i=0;i<10;++i){
                phat.push_back(-one);//error flag
                pNeutral.push_back(false);}
            qhat=phat; qNeutral=pNeutral;
			int numMaxima=0;
            std::string dummyErrorStr("dummyErrorString");
            MaximizePopMeanFitnessPandQv2(wAABB,wAABb,WAAbb,wAaBB,wAaBb,wAabb,waaBB,waaBb,waabb,wBarMax,phat,qhat,pNeutral,qNeutral,numMaxima,dummyErrorStr);
		 
		 	cout<<"wbar="<<wBarMax<<endl;
			for(int m=0;m<numMaxima;++m){
				cout<<"phat="<<phat[m]<<", qhat="<<qhat[m]<<endl;}
		 	cout<<endl;
			return 0;
	
//    DetermineViableMismatchCombinations(3);
*/
	
	int bitstringLen = 4;
    int Ntfsat_int=10;//Ntfsat->{10,30,100,300,1000} ->set the deltaG values too
    int PoptLow=0, PoptSteps=1000, PoptStepSize=5;
    int PoptHigh=1000;
    bool saveAllSolutions=false;//if false, then only the summary table is saved
    bool runningInSegments=false;
    bool runUsingThreads=true;
    bool splitSinglePoptRun=false;
    typeOfModelToRun modelToRun=allSites;
    string outputSummaryFileDesignator, outputSummaryFileHeaderDesignator;
    bool printSeparateHeaderFile=false;//use this for concatenating files outside
//    std::string bb="3";//debugging
//    outputSummaryFileHeaderDesignator="_b"+bb+"_header";
    long double omega= (long double)0.05;//0.05, 0.2, 0.005
	long double NtfsatPerAllele;
	long double deltaG1dosage = (long double) -7.90249;//corresponding deltaG->{-3.89182, -5.66643, -7.90249, -10.0478, -12.4372}
	long double deltaG1 = (long double) -7.86365;//correspondng deltaG->{-3.58352, -5.54518, -7.86365, -10.0346, -12.4332}
    uint64_t maxBitstringVal = uint64_t(pow(2,bitstringLen));
    uint64_t tf0Low=0, tf0High=maxBitstringVal-1;
    
    //debugging splitPoptRun:
//    splitSinglePoptRun=true;
//    tf0Low=2;
//    tf0High=3;
//    PoptLow=PoptHigh=195;
                for(int i=0;i<argc;++i){
                    std::cout<<"argv["<<i<<"]="<<argv[i]<<std::endl;}
    
    if(argc>=8){//get parameters off the input line
            std::string b=argv[1];
            std::string n=argv[2];
            std::string o=argv[3];
            std::string pl=argv[4];
            std::string ph=argv[5];
            std::string ps=argv[6];
            std::string pss=argv[7];
            std::string tf0LowStr, tf0HighStr;
            std::string::size_type sz;
            bitstringLen=std::stoi(b,&sz);
            maxBitstringVal = uint64_t(pow(2,bitstringLen));
            Ntfsat_int=std::stoi(n,&sz);
            omega=std::stold(o,&sz);
            PoptLow=std::stoi(pl,&sz);
            PoptLow = MAX(0,PoptLow);
            PoptHigh=std::stoi(ph,&sz);
            PoptSteps=std::stoi(ps,&sz);
            PoptHigh = MIN(PoptHigh,PoptSteps);
            PoptStepSize=std::stoi(pss,&sz);
            runningInSegments=true;
            saveAllSolutions=false;//should be anyway
            outputSummaryFileHeaderDesignator="_b"+b+"_Ntf"+n+"_header";
            float poptlow = (float)PoptLow/(float)PoptSteps;
            float popthigh = (float)PoptHigh/(float)PoptSteps;
            stringstream sspl; sspl << poptlow;
            outputSummaryFileDesignator+="_b"+b+"_Ntf"+n+"_Popt"+sspl.str();
            if(PoptHigh>PoptLow){
                stringstream ssph; ssph << popthigh;
                outputSummaryFileDesignator+="to"+ssph.str();}
            if(!(PoptLow==0 && PoptHigh==PoptSteps)){//if it's not the full range, assume we're running batch files to concatenate later
                printSeparateHeaderFile=true;}
            if(argc==10){//splitting a single Popt run
                splitSinglePoptRun=true;
                PoptHigh=PoptLow;//in case this is set for >1 Popt value
                tf0LowStr=argv[8];
                tf0HighStr=argv[9];
                std::cout<<"tf0LowStr="<<tf0LowStr<<", tf0HighStr="<<tf0HighStr<<std::endl;
                tf0Low=(uint64_t) stoi(tf0LowStr,&sz);
                tf0High=(uint64_t) stoi(tf0HighStr,&sz);
                std::cout<<"tf0Low="<<tf0Low<<", tf0High="<<tf0High<<std::endl;
                tf0Low=MIN(tf0Low,(uint64_t)maxBitstringVal-1);
                tf0High=MIN(tf0High,(uint64_t)maxBitstringVal-1);
                std::cout<<"tf0LowStr="<<tf0LowStr<<", tf0HighStr="<<tf0HighStr<<std::endl;
                if(tf0High<tf0Low){uint64_t temp=tf0Low; tf0Low=tf0High; tf0High=temp;}//swap them
                outputSummaryFileDesignator+="_tf"+tf0LowStr;
                if(tf0High>tf0Low){
                    outputSummaryFileDesignator+="to"+tf0HighStr;}
                }//argc==10
            }//argc >=8
        else{//use default file designators
            tf0Low=MIN(tf0Low,maxBitstringVal-1);
            tf0High=MIN(tf0High,maxBitstringVal-1);
            std::stringstream bb, nn, oo, tflow, tfhigh, pp;
            bb<<bitstringLen;
            nn<<Ntfsat_int;
            oo<<omega;
            tflow<<tf0Low;
            tfhigh<<tf0High;
            float poptlow = (float)PoptLow/(float)PoptSteps;
            poptlow = ROUND(poptlow,3);
            pp<<poptlow;

            std::string b = bb.str();
            std::string n = nn.str();
            std::string o = oo.str();
            std::string tL = tflow.str();
            std::string tH = tfhigh.str();
            std::string p = pp.str();
            outputSummaryFileHeaderDesignator="_b"+b+"_Ntf"+n+"_header";
            outputSummaryFileDesignator+="_b"+b+"_Ntf"+n+"_o"+o;
            if(splitSinglePoptRun){
                outputSummaryFileDesignator+=+"_Popt"+p+"_tf"+tL;
                if(tf0High>tf0Low){
                    outputSummaryFileDesignator+="to"+tH;}
                }
            }
//	long double deltaG1dosage = (long double) -7.90249;//corresponding deltaG->{-3.89182, -5.66643, -7.90249, -10.0478, -12.4372}
//	long double deltaG1 = (long double) -7.86365;//correspondng deltaG->{-3.58352, -5.54518, -7.86365, -10.0346, -12.4332}
    NtfsatPerAllele=(long double)Ntfsat_int;
    switch(Ntfsat_int){
        case 10:
            deltaG1dosage=-3.89182; deltaG1=-3.58352;
            break;
        case 30:
            deltaG1dosage=-5.66643; deltaG1=-5.54518;
            break;
        case 300:
            deltaG1dosage=-10.0478; deltaG1=-10.0346;
            break;
        case 1000:
            deltaG1dosage=-12.4372; deltaG1=-12.4332;
            break;
        default: //case 100
            deltaG1dosage=-7.90249; deltaG1=-7.86365;
          break;
        }
	SimplestRegPathIndividual minmaxIndiv;
	long double minExpression=zero,maxExpression=zero;
	minmaxIndiv.CalculateMinMaxExpression(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,
			minExpression,maxExpression,modelToRun);
	
//	long double Popt= (long double)0.60;

/*	//debug
	SimplestRegPathIndividual indivAABB,indivAABb,indivAAbb,indivAaBB,indivAaBb,indivAabb,indivaaBB,indivaaBb, indivaabb;
	int mDoseTf0=1,mDoseTf1=2,mTf0cis0=3,mTf0cis1=4,mTf1cis0=5,mTf1cis1=6;
	indivAaBb.SetGenotype(0,0,0,mDoseTf0);indivAaBb.SetGenotype(1,1,0,mDoseTf1);
	indivAaBb.SetGenotype(-1,0,0,mTf0cis0);indivAaBb.SetGenotype(-1,0,1,mTf0cis1);
	indivAaBb.SetGenotype(-1,1,0,mTf1cis0);indivAaBb.SetGenotype(-1,1,1,mTf1cis1);
	
	indivAABB.TFdosage_[0]=indivAABB.TFdosage_[1]=mDoseTf0;
	indivAABB.mTF01cis01_[0][0]=indivAABB.mTF01cis01_[0][1]=mTf0cis0;
	indivAABB.mTF01cis01_[1][0]=indivAABB.mTF01cis01_[1][1]=mTf0cis0;
	
	indivAABb.TFdosage_[0]=indivAABb.TFdosage_[1]=mDoseTf0;
	indivAABb.mTF01cis01_[0][0]=indivAABb.mTF01cis01_[1][0]=mTf0cis0;
	indivAABb.mTF01cis01_[0][1]=indivAABb.mTF01cis01_[1][1]=mTf0cis1;
	
	indivAAbb.TFdosage_[0]=indivAAbb.TFdosage_[1]=mDoseTf0;
	indivAAbb.mTF01cis01_[0][0]=indivAAbb.mTF01cis01_[1][0]=mTf0cis1;
	indivAAbb.mTF01cis01_[0][1]=indivAAbb.mTF01cis01_[1][1]=mTf0cis1;
	
	indivAaBB.TFdosage_[0]=mDoseTf0;
	indivAaBB.TFdosage_[1]=mDoseTf1;
	indivAaBB.mTF01cis01_[0][0]=indivAaBB.mTF01cis01_[0][1]=mTf0cis0;
	indivAaBB.mTF01cis01_[1][0]=indivAaBB.mTF01cis01_[1][1]=mTf1cis0;
	
	indivAabb.TFdosage_[0]=mDoseTf0;
	indivAabb.TFdosage_[1]=mDoseTf1;
	indivAabb.mTF01cis01_[0][0]=indivAabb.mTF01cis01_[0][1]=mTf0cis1;
	indivAabb.mTF01cis01_[1][0]=indivAabb.mTF01cis01_[1][1]=mTf1cis1;
	
	indivaaBB.TFdosage_[0]=indivaaBB.TFdosage_[1]=mDoseTf1;
	indivaaBB.mTF01cis01_[0][0]=indivaaBB.mTF01cis01_[0][1]=mTf1cis0;
	indivaaBB.mTF01cis01_[1][0]=indivaaBB.mTF01cis01_[1][1]=mTf1cis0;
	
	indivaaBb.TFdosage_[0]=indivaaBb.TFdosage_[1]=mDoseTf1;
	indivaaBb.mTF01cis01_[0][0]=indivaaBb.mTF01cis01_[1][0]=mTf1cis0;
	indivaaBb.mTF01cis01_[0][1]=indivaaBb.mTF01cis01_[1][1]=mTf1cis1;
	
	indivaabb.TFdosage_[0]=indivaabb.TFdosage_[1]=mDoseTf1;
	indivaabb.mTF01cis01_[0][0]=indivaabb.mTF01cis01_[1][0]=mTf1cis1;
	indivaabb.mTF01cis01_[0][1]=indivaabb.mTF01cis01_[1][1]=mTf1cis1;
	
	//confirm format in debugger
	std::string mAABB=indivAABB.mismatchStringMathematicaFormat();
	std::string mAABb=indivAABb.mismatchStringMathematicaFormat();
	std::string mAAbb=indivAAbb.mismatchStringMathematicaFormat();
	std::string mAaBB=indivAaBB.mismatchStringMathematicaFormat();
	std::string mAaBb=indivAaBb.mismatchStringMathematicaFormat();
	std::string mAabb=indivAabb.mismatchStringMathematicaFormat();
	std::string maaBB=indivaaBB.mismatchStringMathematicaFormat();
	std::string maaBb=indivaaBb.mismatchStringMathematicaFormat();
	std::string maabb=indivaabb.mismatchStringMathematicaFormat();

*/



//    std::vector<SimplestRegPathIndividual>::iterator i_bestIndivList;
//    std::vector<long double>::iterator i_bestLDph, i_bestLDf;
	std::vector<long double> bestPhenotypes, bestFitnesses, bestMeanPhenotypes;
//	std::vector<std::string>::iterator i_bestStg,i_bestStm,i_bestHt, i_bestLDp, i_bestLDq;
	std::vector<std::string> bestGtypes, bestMismatchGtypes, bestHetType, bestPs, bestQs;
	std::string bestMismatchGtype,pBest,qBest,phBest;
    
    std::string outputfileAllSolutionsName("fitnessOverdominanceDataTable");
    outputfileAllSolutionsName += outputSummaryFileDesignator+".txt";
    char* outputfileAllSolutionsNameStr = new char[outputfileAllSolutionsName.length()+1];
    std::strcpy(outputfileAllSolutionsNameStr,outputfileAllSolutionsName.c_str());
    std::fstream outputfileAllSolutions;
    FitnessMaximaBitstringSolutions wBarMaxAllSolutions;//(splitSinglePoptRun,tf0Low,tf0High);
    if(saveAllSolutions){
        outputfileAllSolutions.open(outputfileAllSolutionsNameStr,std::fstream::out);
        outputfileAllSolutions.close();
        outputfileAllSolutions.open(outputfileAllSolutionsNameStr,std::fstream::app);
        wBarMaxAllSolutions.PrintSolutionTableHeader(outputfileAllSolutions);}

    std::string outputfileSolutionSummariesName("fitnessOverdomSummaryTable");
    outputfileSolutionSummariesName += outputSummaryFileDesignator+".txt";
    char* outputfileSolutionSummariesNameStr = new char[outputfileSolutionSummariesName.length()+1];
    std::strcpy(outputfileSolutionSummariesNameStr,outputfileSolutionSummariesName.c_str());
    std::fstream outputfileSolutionSummaries,outputfileSolutionSummariesHeader;
    outputfileSolutionSummaries.open(outputfileSolutionSummariesNameStr,std::fstream::out);//overwrites
    outputfileSolutionSummaries.close();
    outputfileSolutionSummaries.open(outputfileSolutionSummariesNameStr,std::fstream::app);

    FitnessMaximaSolutionSets summariesOfSolutions(splitSinglePoptRun,tf0Low,tf0High);
    if(printSeparateHeaderFile){
            std::string outputfileSolutionSummariesHeaderName("fitnessOverdomSummaryTable");
            outputfileSolutionSummariesHeaderName += outputSummaryFileHeaderDesignator+".txt";
            char* outputfileSolutionSummariesHeaderNameStr = new char[outputfileSolutionSummariesHeaderName.length()+1];
            std::strcpy(outputfileSolutionSummariesHeaderNameStr,outputfileSolutionSummariesHeaderName.c_str());
            outputfileSolutionSummariesHeader.open(outputfileSolutionSummariesHeaderNameStr,std::fstream::out);//overwrites
            outputfileSolutionSummariesHeader.close();
            outputfileSolutionSummariesHeader.open(outputfileSolutionSummariesHeaderNameStr,std::fstream::app);
            delete [] outputfileSolutionSummariesHeaderNameStr; outputfileSolutionSummariesHeaderNameStr=NULL;
            summariesOfSolutions.PrintHeaderLine(outputfileSolutionSummariesHeader);
            outputfileSolutionSummariesHeader.close();
            }
        else{
            summariesOfSolutions.PrintHeaderLine(outputfileSolutionSummaries);}


    aTime timer;
    for(int i=PoptLow;i<=PoptHigh;i+=PoptStepSize){//+=5 normally
        long double Popt = (long double)i/(long double)PoptSteps;
        switch(modelToRun){
            case dosageOnly:
                MaximizeUsingBitstringsDosageOnly(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression,
                    Popt,omega,wBarMaxAllSolutions,outputfileAllSolutions,saveAllSolutions,
                    summariesOfSolutions,outputfileSolutionSummaries);
                break;
            case tfProductOnly:
                MaximizeUsingBitstringsTFproductOnly(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression,
                    Popt,omega,wBarMaxAllSolutions,outputfileAllSolutions,saveAllSolutions,
                    summariesOfSolutions,outputfileSolutionSummaries);
                break;
            case cisOnly:
                MaximizeUsingBitstringsCisOnly(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression,
                    Popt,omega,wBarMaxAllSolutions,outputfileAllSolutions,saveAllSolutions,
                    summariesOfSolutions,outputfileSolutionSummaries);
                break;
            case tfOnly: //dosage & product; always uses threads
                MaximizeUsingBitstringsTFOnly(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression,
                    Popt,omega,wBarMaxAllSolutions,outputfileAllSolutions,saveAllSolutions,
                    summariesOfSolutions,outputfileSolutionSummaries);
                break;
            default: //allSites
                if(runUsingThreads){
                        MaximizeUsingBitstringsThreadableAllCis(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,
                            minExpression,maxExpression,Popt,omega,wBarMaxAllSolutions,outputfileAllSolutions,saveAllSolutions,
                            summariesOfSolutions,outputfileSolutionSummaries,splitSinglePoptRun,tf0Low,tf0High);
                        }
                    else{
                        MaximizeUsingBitstrings(bitstringLen,NtfsatPerAllele,deltaG1dosage,deltaG1,minExpression,maxExpression,
                            Popt,omega,wBarMaxAllSolutions,outputfileAllSolutions,saveAllSolutions,summariesOfSolutions,
                                outputfileSolutionSummaries,splitSinglePoptRun,tf0Low,tf0High);}
            }//switch modelToRun


        std::string *elapsed=timer.HMS_elapsed();
        coutLock.lock(); std::cout<<"elapsed = "<<*elapsed<<std::endl; coutLock.unlock();
        delete elapsed;elapsed=NULL;
        }
    outputfileAllSolutions.close();
    outputfileSolutionSummaries.close();
    delete [] outputfileAllSolutionsNameStr;
    delete [] outputfileSolutionSummariesNameStr;
    outputfileAllSolutionsNameStr=outputfileSolutionSummariesNameStr=NULL;
    return 0;
}//main

