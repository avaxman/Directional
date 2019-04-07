#ifndef DIRECTIONAL_RANGEHELPERS_H
#define DIRECTIONAL_RANGEHELPERS_H
#include <vector>
namespace Helpers{
    template<typename T>
    inline std::vector<T> Constant(int size, T t){
        std::vector<T> els;
		els.reserve(size);
        for(int i =0 ; i < size; i++) els.push_back(t);
		return els;
    }
	inline void intRange(int start, int endExcl, std::vector<int>& target)
	{
		target.resize(endExcl - start);
		for (int i = start; i < endExcl; i++) target[i - start] = i;
	}

	void assignRange(int start, int endExcl, std::vector<int>& target)
	{
		target.resize(endExcl - start);
		for (int i = 0; i < endExcl - start; i++)target[i] = start + i;
	}
	void assignWrappedRange(int start, int endExcl, int wrap, std::vector<int>& target)
    {
		target.resize(endExcl - start);
	    for(int v = start; v < endExcl; v++)
	    {
			int val = v;
			if (val < 0) val += wrap;
			if (val > wrap) val -= wrap;
			target[v-start] = val;
	    }
    }
	void assignDecreasingRange(int start, int endExcl, std::vector<int>& target)
	{
		target.resize(endExcl - start);
		for (int v = endExcl - 1, i = 0; v >= start; v--, i++)target[i] = v;
	}
	void assign(std::vector<double>& target, const std::vector<double>& values, double factor = 1.0)
    {
		target = std::vector<double>(values.size());
		for (int i = 0; i < values.size(); i++) target[i] = values[i] * factor;
    }
}
#endif