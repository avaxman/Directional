//
// Created by Amir Vaxman on 20.04.24.
//

#ifndef DIRECTIONAL_ENUMBER_DEFINITIONS_H
#define DIRECTIONAL_ENUMBER_DEFINITIONS_H

#include <set>
#include <math.h>
#include <vector>
#include <queue>
#include <cassert>
#include <algorithm>
#include <utility>
#include <Eigen/Sparse>
#include <Eigen/Dense>
//#include <BigInt.hpp>
//#include <functions/math.hpp>
//#include <operators/arithmetic_assignment.hpp>
#include <directional/BigInteger.h>


namespace directional{

typedef BigInteger EInt;

class ENumber{
public:
  EInt num, den;
  bool simple;  //whether number is simplified
  
  ENumber(){num=0; den=1;}
  
  /*ENumber(const double number, const double resolution=10e-9){
   simple=true;
   }*/
  
  ~ENumber(){}
  
  //Computed a continuous fraction approximation to a given tolerance>0
  ENumber(double x, double tol) {
    long long prevNum(0), prevDen(1), currNum(1), currDen(0);
    double fraction = x;
    
    
    while (true) {
      long long integerPart = static_cast<long long>(fraction);
      double remainder = fraction - integerPart;
      
      // Update current numerator and denominator
      long long newNum = integerPart * currNum + prevNum;
      long long newDen = integerPart * currDen + prevDen;
      
      double approximation = static_cast<double>(newNum) / static_cast<double>(newDen);
      
      if (std::abs(approximation - x) <= tol) {
        if (newDen < 0) {newDen=-newDen; newNum = -newNum;}
        num = BigInteger(newNum);
        den = BigInteger(newDen);
        break;
      }
      
      // Update previous and current numerators/denominators
      prevNum = currNum;
      prevDen = currDen;
      currNum = newNum;
      currDen = newDen;
      
      if (remainder == 0) {
        break;
      }
      
      fraction = 1.0 / remainder;
    }
    //std::cout<<"x: "<<x<<std::endl;
    //std::cout<<"approximation: "<<to_double()<<std::endl;
    //std::cout<<"num: "<<num.to_string()<<std::endl;
    //std::cout<<"den: "<<den.to_string()<<std::endl;
  }
  
  ENumber(const EInt _num, const EInt _den, const bool toSimplify=true):num(_num), den(_den), simple(true){
    assert("ENumber(): denominator is zero!" && den!=0);
    if (toSimplify) simplify();
    if (den<0){
      den = -den;
      num = -num;}
  }
  
  ENumber(const EInt _num){
    num = _num;
    den = 1;
  }
  
  void simplify(){
    if (num==0) {den=1; return;}
    EInt common = gcd(num, den);
    num/=common;
    den/=common;
    if (den<0){
      den = -den;
      num = -num;
    }
    
    simple=true;
  }
  
  ENumber operator=(const ENumber& e2){
    num = e2.num;
    den = e2.den;
    return *this;
  }
  
  ENumber operator+(const ENumber& b2) const{
    if (this->den == b2.den)
      return ENumber(this->num + b2.num, this->den);
    else
      return ENumber((this->den * b2.num + this->num * b2.den), (this->den * b2.den));
  }
  
  
  ENumber operator+=(const ENumber& b2){
    *this = *this + b2;
    return *this;
  }
  
  ENumber operator-(const ENumber& b2) const{
    if (this->den == b2.den)
      return ENumber(this->num - b2.num, this->den);
    else
      return ENumber((this->num * b2.den - this->den * b2.num), (this->den * b2.den));
  }
  
  ENumber operator-() const{
    return ENumber(-num, den, false);
  }
  
  ENumber operator*(const ENumber& b2) const{
    ENumber cross1(this->num, b2.den);
    ENumber cross2(b2.num, this->den);
    return ENumber(cross1.num * cross2.num, cross1.den * cross2.den, false);
  }
  
  ENumber operator/(const ENumber& b2) const{
    assert("ENumber division by zero!" && b2.num!=0);
    ENumber cross1(this->num, b2.num);
    ENumber cross2(b2.den, this->den);
    return ENumber(cross1.num * cross2.num, cross1.den * cross2.den, false);
  }
  
  ENumber operator/=(const ENumber& b2){
    *this = *this / b2;
    return *this;
  }
  
  //Numbers are simple now because they are constructed to be simple
  bool operator==(const ENumber& b2) const{
    //return (this->num*b2.den == this->den*b2.num);
    return ((this->num == b2.num)&&(this->den==b2.den));
  }
  
  bool operator!=(const ENumber& b2) const{
    //return (this->num*b2.den != this->den*b2.num);
    return ((this->num != b2.num)||(this->den!=b2.den));
  }
  
  bool operator>=(const ENumber& b2) const{
    if ((b2.num<=0)&&(this->num>=0))
      return true;
    if ((b2.num>0)&&(this->num<0))
      return false;
    return (this->num*b2.den >= this->den*b2.num);
  }
  
  bool operator<=(const ENumber& b2) const{
    if ((b2.num>=0)&&(this->num<=0))
      return true;
    if ((b2.num<0)&&(this->num>0))
      return false;
    return (this->num*b2.den <= this->den*b2.num);
  }
  
  bool operator>(const ENumber& b2) const{
    if ((b2.num<0)&&(this->num>0))
      return true;
    if ((b2.num>0)&&(this->num<0))
      return false;
    return (this->num*b2.den > this->den*b2.num);
  }
  
  bool operator<(const ENumber& b2) const{
    if ((b2.num>0)&&(this->num<0))
      return true;
    if ((b2.num<0)&&(this->num>0))
      return false;
    return (this->num*b2.den < this->den*b2.num);
  }
  
  ENumber abs() const{
    return ENumber((num >0 ? num : -num), den, false);
  }
  
  //getting one by one digits
  long double to_double(const int maxDigits = 12) const{
    EInt currNum = num;
    EInt currDen = den;
    std::string mantissa = (num > 0 ? "+" : "-");
    currNum = (num > 0 ? num : -num);
    EInt quotient = currNum/currDen;
    EInt currRem = currNum - quotient * currDen;
    mantissa+=quotient.to_string() + ".";
    for (int i=0;i<maxDigits;i++)
    {
      currNum = currRem * 10;
      quotient = currNum/currDen;
      currRem = currNum - quotient * currDen;
      mantissa += quotient.to_string();
    }
    //std::ostringstream oss;
    //oss << std::showpos << 0;  // std::showpos forces a + sign for positive numbers
    //std::string exponent = oss.str();
    double result = std::stod(mantissa);
    return result;
  }
  
  /*bool operator==(const BigInt& n2) const{
   return ((this->num == b2.num)&&(this->den==b2.den));
   }
   
   bool operator!=(const ENumber& b2) const{
   return ((this->num != b2.num)||(this->den!=b2.den));
   }*/
};

template <size_t Size>
class EVector {
public:
  EVector() : data(Size) {}
  
  // Other methods can be added as needed
  const ENumber& operator[](const size_t index) const {
    return data[index];
  }
  
  ENumber& operator[](const size_t index) {
    return data[index];
  }
  
  EVector<Size> operator+(const EVector<Size>& ev) const{
    EVector<Size> newVec;
    for (int i=0;i<Size;i++)
      newVec.data[i]=data[i]+ev.data[i];
    
    return newVec;
  }
  EVector<Size> operator-() const{
    EVector<Size> newVec;
    for (int i=0;i<Size;i++)
      newVec.data[i]=-data[i];
    
    return newVec;
  }
  EVector<Size> operator-(const EVector<Size>& ev) const{
    EVector<Size> newVec;
    for (int i=0;i<Size;i++)
      newVec.data[i]=data[i]-ev.data[i];
    
    return newVec;
  }
  
  EVector<Size> operator*(const ENumber s) const{
    EVector<Size> newVec;
    for (int i=0;i<Size;i++)
      newVec.data[i]=data[i]*s;
    
    return newVec;
  }
  
  bool operator==(const EVector<Size>& ev) const{
    bool equal=true;
    for (int i=0;i<Size;i++)
      equal = equal & (ev.data[i]==data[i]);
    return equal;
  }
  
  //for the sake of sorting
  bool operator<(const EVector<Size>& ev) const{
    for (int i=0;i<Size;i++)
      if (data[i]!=ev.data[i])
        return data[i]<ev.data[i];
    return false;
  }
  
  EVector(const std::initializer_list<ENumber>& args){
    data.insert(data.end(), args.begin(), args.end());
  }
  
  EVector operator=(const EVector<Size>& evec){
    data=evec.data;
    return *this;
  }
  
  template<size_t _Size>
  friend std::ostream& operator<<(std::ostream& os, const EVector<_Size>& evec);
  
  /*void canonicalize(){
   for (int i=0;i<Size;i++)
   data[i].canonicalize();
   }*/
  
  //protected:
  std::vector<ENumber> data;
};

template<size_t Size>
EVector<Size> operator*(ENumber scalar, const EVector<Size>& vec) {
  return vec * scalar; // Leverage the previous operator*
}

template<size_t Size>
std::ostream& operator<<(std::ostream& os, const EVector<Size>& evec) {
  os << "(";
  for (int i=0;i<Size-1;i++)
    os<< evec[i].to_double()<<",";
  os<<evec[Size-1].to_double()<<")";
  return os;
}



typedef EVector<2> EVector2;
typedef EVector<3> EVector3;

struct Segment2 {
public:
  EVector2 source, target;
  Segment2(const EVector2& _source, const EVector2& _target){
    source=_source; target=_target;
  }
  
  Segment2 operator=(const Segment2& seg2){
    source = seg2.source;
    target = seg2.target;
    return *this;
  }
  
  Segment2(){}
  
  friend std::ostream& operator<<(std::ostream& os, const Segment2& seg);
};

std::ostream& operator<<(std::ostream& os, const Segment2& seg) {
  os << "Segment2(" << seg.source << "->" << seg.target << ")";
  return os;
}

struct Line2{
public:
  EVector2 point, direction;
  Line2(const EVector2& _point, const EVector2& _direction){
    point=_point; direction=_direction;
  }
  
  friend std::ostream& operator<<(std::ostream& os, const Line2& seg);
};

std::ostream& operator<<(std::ostream& os, const Line2& line) {
  os << "Line2(" << line.point << " + " << line.direction << ")";
  return os;
}


struct LinePencil{
  int minIsoValue, maxIsoValue;
  EVector2 direction;
  EVector2 p0, pVec;  //p0 is the position of the line at minValue. pVec is the vector between the origins of the lines
};


ENumber squaredDistance(const EVector3& v1, const EVector3& v2){
  ENumber sd(0);
  //assert("vectors are of different dimensions!" && v1.size()==v2.size());
  for (int i=0;i<3;i++)
    sd+=(v1[i]-v2[i])*(v1[i]-v2[i]);  //maybe it's not efficient
  
  return sd;
}

//produces y = M*x
void exactSparseMult(const Eigen::SparseMatrix<int> M, const std::vector<ENumber>& x,std::vector<ENumber>& y){
  y.resize(M.rows());
  
  for (int i=0;i<y.size();i++)
    y[i]=ENumber(0);
  
  for (int k=0; k<M.outerSize(); ++k)
    for (Eigen::SparseMatrix<int>::InnerIterator it(M,k); it; ++it){
      /*if (it.row()==219){
       std::cout<<"it.value(): "<<it.value()<<std::endl;
       std::cout<<"x[it.col()]: "<<x[it.col()].to_double()<<std::endl;
       std::cout<<"ENumber((long)it.value())*x[it.col()]: "<<(ENumber((long)it.value())*x[it.col()]).to_double()<<std::endl;
       std::cout<<"y[it.row()]: "<<y[it.row()].to_double()<<std::endl;
       }*/
      y[it.row()]+=ENumber((long)it.value())*x[it.col()];
      
    }
}

void exactDenseMult(const Eigen::MatrixXi &nM, const Eigen::MatrixXi& dM, const std::vector<ENumber>& x, std::vector<ENumber>& y)
{
  y.resize(nM.rows());
  for (int i=0;i<y.size();i++)
    y[i]=ENumber(0);
  for (int i=0;i<nM.rows();i++)
    for (int j=0;j<nM.cols();j++)
      y[i]+=x[j]*ENumber(nM(i,j), dM(i,j));
}

//This assumes components is already resized to the correct |v|
//not very efficient but probably not terrible
int connectedComponents(const std::vector<std::pair<int,int>>& matches, std::vector<int>& components)
{
  for (int i=0;i<components.size();i++)
    components[i]=-1;
  
  std::vector<std::vector<int>> VV(components.size());
  for (int i=0;i<matches.size();i++){
    VV[matches[i].first].push_back(matches[i].second);
    VV[matches[i].second].push_back(matches[i].first);
  }
  
  std::deque<int> nextVertexQueue;
  for (int i=0;i<components.size();i++)
    nextVertexQueue.push_front(i);
  
  int numComponents=0;
  while (!nextVertexQueue.empty()) {
    int nextVertex=nextVertexQueue.front();
    nextVertexQueue.pop_front();
    if (components[nextVertex]==-1) {  //first components
      std::cout << "New component " << numComponents << " seed vertex " << nextVertex << std::endl;
      components[nextVertex] = numComponents++;
    }
    
    //Otherwise, doing DFS on edges
    for (int i=0;i<VV[nextVertex].size();i++){
      if (components[VV[nextVertex][i]]==-1){
        components[VV[nextVertex][i]]=components[nextVertex];
        std::cout<<"adding vertex "<<VV[nextVertex][i]<<" to component "<<components[nextVertex]<<std::endl;
        nextVertexQueue.push_front(VV[nextVertex][i]);
      } else assert(components[VV[nextVertex][i]]==components[nextVertex]);
    }
  }
  return numComponents;
}

int line_line_intersection(const Line2& line1,
                           const Line2& line2,
                           ENumber& t1,
                           ENumber& t2){
  //std::cout<<"Computing intersection between line: "<<line1<<" and line "<<line2<<std::endl;
  ENumber denom = (line1.direction[0]*line2.direction[1]-line1.direction[1]*line2.direction[0]);
  if (denom==ENumber(0)){
    EVector2 pointVec = line1.point-line2.point;
    //std::cout<<"lines are parallel"<<std::endl;
    return  (pointVec[0]*line1.direction[1]-pointVec[1]*line1.direction[0]==ENumber(0) ? 2 : 0);
  }
  
  //std::cout<<"exact line2.point[0]-line1.point[0]: "<<(line2.point[0]-line1.point[0]).to_double()<<std::endl;
  //std::cout<<"double line2.point[0]-line1.point[0]: "<<(line2.point[0].to_double()-line1.point[0].to_double())<<std::endl;
  
  //std::cout<<"exact (line2.point[0]-line1.point[0])*(line2.direction[1]): "<<((line2.point[0]-line1.point[0])*(line2.direction[1])).to_double()<<std::endl;
  //std::cout<<"double (line2.point[0]-line1.point[0])*(line2.direction[1]):  "<<(line2.point[0].to_double()-line1.point[0].to_double())*(line2.direction[1].to_double())<<std::endl;
  
  t1 = ((line2.point[0]-line1.point[0])*(line2.direction[1])-(line2.point[1]-line1.point[1])*(line2.direction[0]))/denom;
  t2 = ((line2.point[0]-line1.point[0])*(line1.direction[1])-(line2.point[1]-line1.point[1])*(line1.direction[0]))/denom;
  //std::cout<<"t1, t2: "<<t1.to_double()<<","<<t2.to_double()<<std::endl;
  //std::cout<<"line1.point+t1*line1.direction: "<<line1.point+t1*line1.direction<<std::endl;
  //std::cout<<"line2.point+t2*line2.direction: "<<line2.point+t2*line2.direction<<std::endl;
  EVector2 diff = (line1.point+t1*line1.direction) - (line2.point+t2*line2.direction);
  //diff.canonicalize();
  assert("line_line_intersection is wrong!" && diff == EVector2());
  //std::cout<<"lines intersect at "<<(line1.point+t1*line1.direction)<<std::endl;
  //std::cout<<"parameters: t1:"<<t1.get_d()<<", t2: "<<t2.get_d()<<std::endl;
  return 1;
  
}


//returns a generator for the grid of intersections, parameterized by p00 + pVec1*isoValue1 + pVec2*isoValue2,
//txp00 is the t(1 or 2) of the p00 point in each respective line
//dtx (1 or 2) is the dt going along each pVecx (1 or 2)
//for iso1 and iso2 in the respective line pencil ranges
//result = 2 is only acceptable if |lp2| = 1, not handling parallel full line pencils (shouldn't be unless the parameterization is degenerate).
inline int linepencil_intersection(const LinePencil& lp1,
                                   const LinePencil& lp2,
                                   directional::EVector2& p00,
                                   directional::EVector2& pVec1,
                                   directional::EVector2& pVec2,
                                   ENumber& t1p00,
                                   ENumber& t2p00,
                                   ENumber& dt1,
                                   ENumber& dt2,
                                   EInt& iso1Overlap){   //overlap line in lp1 matching that of lp2. ignored unless result==2
  
  
  ENumber denom = (lp1.direction[0]*lp2.direction[1]-lp1.direction[1]*lp2.direction[0]);
  if (denom==ENumber(0)){  //the two pencils are parallel; looking which line overlaps lp2
    assert("linepencil_intersection: parallel line pencils must have |lp2|=1; parameterization must be degenerate! " && lp2.isoMinValue==lp2.isoMaxValue);
    EVector2 pointVec = lp2.p0-lp1.p0;
    ENumber pVecCrossDirection1 = lp1.pVec[0]*lp1.direction[1] - lp1.pVec[1]*lp1.direction[0];
    assert("pVecCrossDirection1 shouldn't be zero! " && pVecCrossDirection1!=ENumber(0));
    ENumber overlapIso = (pointVec[0]*lp1.direction[1]-pointVec[1]*lp1.direction[0]) /pVecCrossDirection1;
    if (overlapIso.den==1){
      iso1Overlap = overlapIso.num;
      return ((iso1Overlap >= lp1.minIsoValue && iso1Overlap <= lp1.maxIsoValue) ? 2 : 0);  //only an overlap if integer and within range of lp1
    } else return 0;  //not an integer; no overlap for certain
  }
  
  //line pencils are next not overlapping; finding grid parameters
  
  //p00 is the intersection of both min isovalues
  t1p00 = ((lp2.p0[0]-lp1.p0[0])*(lp2.direction[1])-(lp2.p0[1]-lp1.p0[1])*(lp2.direction[0]))/denom;
  t2p00 = ((lp2.p0[0]-lp1.p0[0])*(lp1.direction[1])-(lp2.p0[1]-lp1.p0[1])*(lp1.direction[0]))/denom;
  p00 = (lp1.p0+t1p00*lp1.direction);
  //TODO: continue from here.
  //std::cout<<"t1, t2: "<<t1.to_double()<<","<<t2.to_double()<<std::endl;
  //std::cout<<"line1.point+t1*line1.direction: "<<line1.point+t1*line1.direction<<std::endl;
  //std::cout<<"line2.point+t2*line2.direction: "<<line2.point+t2*line2.direction<<std::endl;
  EVector2 diff = (lp1.p0 + t1p00 * lp1.direction) - (lp2.p0 + t2p00 * lp2.direction);
  //diff.canonicalize();
  assert("line_line_intersection is wrong!" && diff == EVector2());
  //std::cout<<"lines intersect at "<<(line1.point+t1*line1.direction)<<std::endl;
  //std::cout<<"parameters: t1:"<<t1.get_d()<<", t2: "<<t2.get_d()<<std::endl;
  return 1;
  
}

std::vector<std::pair<ENumber, ENumber>> segment_segment_intersection(const Segment2& seg1,
                                                                      const Segment2& seg2){
  
  ENumber t1, t2;
  //std::cout<<"Computing intersection of "<<seg1<<" and "<<seg2<<std::endl;
  int result = line_line_intersection(Line2(seg1.source, seg1.target-seg1.source),
                                      Line2(seg2.source, seg2.target-seg2.source),t1, t2);
  
  if (result==0) {
    //std::cout<<"supporting lines don't intersect"<<std::endl;
    return std::vector<std::pair<ENumber, ENumber>>(); //no intersection
  }
  
  if (result==1) {  //a single intersection at most; should check t1 and t2
    //std::cout<<"single intersection"<<std::endl;
    if ((t1>=ENumber(0))&&(t1<=ENumber(1))&&(t2>=ENumber(0))&&(t2<=ENumber(1))){
      std::vector<std::pair<ENumber, ENumber>> point(1);
      point[0]=std::pair<ENumber, ENumber>(t1,t2);
      //std::cout<<"Intersecting at parameters "<<point[0].first.get_d()<<","<<point[0].second.get_d()<<std::endl;
      return point;
    } else{
      //std::cout<<"Intersecting out of parameter bounds"<<std::endl;
      return std::vector<std::pair<ENumber, ENumber>>(); //no intersection
    }
  }
  
  if (result==2){  //lines overlap; should check the segments overlap and then return both overlap points (order not important)
    //std::cout<<"Supporting lines overlap"<<std::endl;
    EVector2 vec = seg1.target-seg1.source;
    int axis = (vec[0]!=ENumber(0) ? 0 : 1);
    Segment2 sortSeg1, sortSeg2;
    if (seg1.source[axis]<seg1.target[axis]) sortSeg1=seg1; else sortSeg1=Segment2(seg1.target,seg1.source);
    if (seg2.source[axis]<seg2.target[axis]) sortSeg2=seg2; else sortSeg2=Segment2(seg2.target,seg2.source);
    EVector2 startPoint = (sortSeg1.source[axis]>sortSeg2.source[axis] ? sortSeg1.source : sortSeg2.source);
    EVector2 endPoint = (sortSeg1.target[axis]>sortSeg2.target[axis] ? sortSeg2.target : sortSeg1.target);
    
    //TODO: what is the policy with segment-vertex tangency?
    if (startPoint[axis]<endPoint[axis]) { //there is a (non-zero) intersection
      ENumber startAtSeg1 = (startPoint[axis]-seg1.source[axis])/(seg1.target[axis]-seg1.source[axis]);
      ENumber startAtSeg2 = (startPoint[axis]-seg2.source[axis])/(seg2.target[axis]-seg2.source[axis]);
      ENumber endAtSeg1 = (endPoint[axis]-seg1.source[axis])/(seg1.target[axis]-seg1.source[axis]);
      ENumber endAtSeg2 = (endPoint[axis]-seg2.source[axis])/(seg2.target[axis]-seg2.source[axis]);
      std::vector<std::pair<ENumber, ENumber>> points(2);
      points[0] = std::pair<ENumber, ENumber>(startAtSeg1, startAtSeg2);
      points[1] = std::pair<ENumber, ENumber>(endAtSeg1, endAtSeg2);
      //std::cout<<"Intersecting at parameters "<<points[0].first.to_double()<<","<<points[0].second.to_double()<<" and "<<points[1].first.to_double()<<","<<points[1].second.to_double()<<std::endl;
      return points;
    } else{
      //std::cout<<"parameters don't overlap"<<std::endl;
      return std::vector<std::pair<ENumber, ENumber>>(); //no intersection
    }
    
  }
  
  
}



std::vector<ENumber> line_segment_intersection(const Line2& line,
                                               const Segment2& segment){
  
  //std::cout<<"Computing intersection between line :"<<line<<" and segment "<<segment<<std::endl;
  Line2 segLine(segment.source, segment.target-segment.source);
  ENumber t1, t2;
  int intersectType=line_line_intersection(line, segLine, t1, t2);
  if (intersectType==0) {  //no intersection
    //std::cout<<"No intersection (lines parallel)"<<std::endl;
    return std::vector<ENumber>();
  } else if (intersectType==2) { //the entire segment is contained in the line
    std::vector<ENumber> result(2); result[0]=ENumber(0); result[1]=ENumber(1);
    std::cout<<"Entire segment is contained in line"<<std::endl;
    return result;
  } else { //(intersectType==1)
    if ((t2>=ENumber(0)) && (t2<=ENumber(1))){
      std::vector<ENumber> result(1); result[0]=t1;
      //std::cout<<"Intersecting at line parameter t1:"<<t1.get_d()<<std::endl;
      return result;
    } else {
      //std::cout<<"No intersection  (out of segment parameter)"<<std::endl;
      return std::vector<ENumber>();
    }
  }
}

void line_triangle_intersection(const Line2& line,
                                const std::vector<EVector2> triangle,
                                bool& intEdge,
                                bool& intFace,
                                ENumber& inParam,
                                ENumber& outParam){
  
  inParam = ENumber(3276700);
  outParam = ENumber(-3276700);
  intFace=intEdge=false;
  for (int i=0;i<3;i++){
    Segment2 edgeSegment(triangle[i], triangle[(i+1)%3]);
    std::vector<ENumber> result = line_segment_intersection(line, edgeSegment);
    for (int j=0;j<result.size();j++){
      inParam = (inParam < result[j] ? inParam : result[j]);
      outParam = (outParam > result[j] ? outParam : result[j]);
    }
    if (result.size()==2){
      std::cout<<"line triangle overlap!!!"<<std::endl;
      intEdge=true;
      intFace=false;
      return;
    }
    if (result.size()==1)
      intFace=true;
  }
  if (inParam==outParam)  //intersecting the triangle only by a vertex; ignored
    intFace=intEdge=false;
  /*if (intFace){
   std::cout<<"Intersecting within the face with parameters "<<inParam.get_d()<<"->"<<outParam.get_d()<<std::endl;
   }*/
  /*else if (intEdge){
   std::cout<<"Intersecting an edge with parameters "<<inParam.get_d()<<"->"<<outParam.get_d()<<std::endl;
   } else std::cout<<"Line doesn't intersect triangle"<<std::endl;*/
  
}

void linepencil_triangle_intersection(const LinePencil& lp,
                                      const std::vector<EVector2> triangle,
                                      std::vector<bool>& intEdges,
                                      std::vector<bool>& intFaces,
                                      std::vector<ENumber>& inParams,
                                      std::vector<ENumber>& outParams){
  
  inParams.resize(lp.maxIsoValue-lp.minIsoValue+1);
  outParams.resize(lp.maxIsoValue-lp.minIsoValue+1);
  intEdges.resize(lp.maxIsoValue-lp.minIsoValue+1);
  intFaces.resize(lp.maxIsoValue-lp.minIsoValue+1);
  for (int i=0;i<inParams.size();i++){
    inParams[i] = ENumber(3276700);
    outParams[i] = ENumber(-3276700);
    intEdges[i]=intFaces[i]=false;
  }
  
  directional::EVector2 p00, pVec1, pVec2;
  ENumber t1p00, t2p00, dt1, dt2;
  for (int i=0;i<3;i++){
    LinePencil triEdgePencil;
    triEdgePencil.minIsoValue = triEdgePencil.maxIsoValue = 0;
    triEdgePencil.p0 = triangle[i];
    triEdgePencil.direction = triangle[(i+1)%3] - triangle[i];
    linepencil_intersection(lp, triEdgePencil, p00, pVec1, pVec2, t1p00, t2p00, dt1, dt2);
    
    Segment2 edgeSegment(triangle[i], triangle[(i+1)%3]);
    std::vector<ENumber> result = line_segment_intersection(line, edgeSegment);
    for (int j=0;j<result.size();j++){
      inParam = (inParam < result[j] ? inParam : result[j]);
      outParam = (outParam > result[j] ? outParam : result[j]);
    }
    if (result.size()==2){
      std::cout<<"line triangle overlap!!!"<<std::endl;
      intEdge=true;
      intFace=false;
      return;
    }
    if (result.size()==1)
      intFace=true;
  }
  if (inParam==outParam)  //intersecting the triangle only by a vertex; ignored
    intFace=intEdge=false;
  /*if (intFace){
   std::cout<<"Intersecting within the face with parameters "<<inParam.get_d()<<"->"<<outParam.get_d()<<std::endl;
   }*/
  /*else if (intEdge){
   std::cout<<"Intersecting an edge with parameters "<<inParam.get_d()<<"->"<<outParam.get_d()<<std::endl;
   } else std::cout<<"Line doesn't intersect triangle"<<std::endl;*/
  
}

//according to this: https://math.stackexchange.com/questions/1450498/rational-ordering-of-vectors
ENumber slope_function(const EVector2& vec){
  //predicates might be expensive, so precomputing
  bool x0 = vec[0]>ENumber(0);
  bool y0 = vec[1]>ENumber(0);
  //bool xy = (y0 && vec[1]>vec[0])||(!y0 && vec[1]<=vec[0]);
  bool xy = vec[1].abs() > vec[0].abs();
  
  //std::cout<<"vec: "<<vec<<std::endl;
  
  if (xy){
    if (y0) return (vec[1]-vec[0])/(vec[1]); // case 1
    else return (vec[1]-vec[0])/(vec[1])+ENumber(4);  //case 3
  } else {
    if (x0){
      if (y0) return (vec[1]-vec[0])/(vec[0]);  //case 0
      else return (vec[1]-vec[0])/(vec[0])+ENumber(8); //case 4
    }else{
      return (vec[1]-vec[0])/(vec[0])+ENumber(4); //case 2
    }
  }
}


ENumber signed_face_area(const std::vector<EVector2>& faceVectors)
{
  EVector2 currVertex; currVertex[0]=currVertex[1]=ENumber(0);
  ENumber sfa(0);
  for (int i=0;i<faceVectors.size();i++){
    EVector2 nextVector = faceVectors[i];
    EVector2 nextVertex = currVertex+nextVector;
    sfa = sfa + currVertex[0]*nextVertex[1]-currVertex[1]*nextVertex[0];
    //std::cout<<"curr sfa: "<<sfa.get_d()<<std::endl;
    currVertex=nextVertex;
  }
  return sfa;
}

ENumber triangle_area(const std::vector<EVector2>& tri){
  EVector e12 = tri[1]-tri[0];
  EVector e13 = tri[2]-tri[0];
  return (e12[0]*e13[1]-e13[0]*e12[1])/ENumber(2);
}

void div_mod(const EInt a, const EInt b, EInt& q, EInt& r){
  //mpz_tdiv_qr(q.get_mpz_t(),r.get_mpz_t(),a,b);
  q = a / b;
  r = a - b * q;
}
}







#endif //DIRECTIONAL_GMP_DEFINITIONS_H
