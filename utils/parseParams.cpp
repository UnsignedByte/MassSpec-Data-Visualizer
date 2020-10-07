/*
* @Author: UnsignedByte
* @Date:   03:04:47, 05-Aug-2020
* @Last Modified by:   UnsignedByte
* @Last Modified time: 18:48:12, 06-Oct-2020
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <codecvt>
#include <locale>
#include <regex>

// [[Rcpp::plugins("cpp17")]]

#if __has_include(<mex.hpp>)
#include <mex.hpp>
#include <mexAdapter.hpp>
matlab::data::ArrayFactory factory;

using Value = matlab::data::Array;
using Array = matlab::data::CellArray;
using NamedArray = matlab::data::Array;

namespace helper {
	template <typename Numeric>
	Value createNumeric(Numeric i) {
		return factory.createScalar(i);
	}

	template <typename Str>
	Value createString(Str s) {
		return factory.createCharArray(s);
	}

	Array createList(size_t l) {
		return factory.createCellArray({1, l});
	}

	NamedArray createNamed(std::vector<std::string> names, std::vector<Value> values) {
		NamedArray n = factory.createStructArray({1,1}, names);
		for (int i = 0; i < names.size(); i++){
			n[0][names[i]] = values[i];
		}
		return n;
	}
};
#elif __has_include("Rcpp.h")
#include "Rcpp.h"
using Value = SEXP;
using Array = Rcpp::List;
using NamedArray = Rcpp::List;

namespace helper {
	template <typename Numeric>
	Value createNumeric(Numeric i) {
		return Rcpp::NumericVector::create(i);
	}

	template <typename Str>
	Value createString(Str s) {
		return Rcpp::CharacterVector::create(s);
	}

	Array createList(size_t l) {
		return Rcpp::List(l);
	}

	NamedArray createNamed(std::vector<std::string> names, std::vector<Value> values) {
		Rcpp::List n;
		for (int i = 0; i < names.size(); i++){
			n[names[i]] = values[i];
		}
		return n;
	}
};
#else
// no file, won't compile.
#endif

// trim from start (in place)
static inline void ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
		return !std::isspace(ch);
	}));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
		return !std::isspace(ch);
	}).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
	ltrim(s);
	rtrim(s);
}

// trim from both ends (copying)
static inline std::string trim_copy(std::string s) {
	trim(s);
	return s;
}

struct var {
	int type;
	std::string name;
	double numeric;
	std::string str;
	std::vector<double> numericArr;
	std::vector<std::string> strArr;
};

//Hash a cstr
constexpr unsigned int hash(const char *s, int off = 0) {
	return !s[off] ? 5381 : (hash(s, off+1)*33) ^ s[off];
} 

NamedArray parse(const std::string& fname){
	std::size_t splitter = fname.find_last_of('.');

	std::string type = fname.substr(splitter+1);
	std::string name = fname.substr(0, splitter-1);

	std::vector<std::string> names;
	std::vector<Value> mvalues;

	std::ifstream fin;
	fin.open("params.p");
	if (!fin.good()) {
		std::cout << "Missing params.p file. Created in home directory." << std::endl;
		std::ofstream outfile ("params.p");
		outfile << "# Find params.p formatting in README.md" << std::endl;
		outfile.close();
		fin.open("params.p");
	}

	bool p = 0;

	std::string l;
	while(fin.peek()!=EOF){
		std::getline(fin, l);

		if (l.length()==0 || l.at(0)=='#') continue;

		if (l.at(0) == '[' && l.at(l.length()-1) == ']') {
			if (l.compare(1, l.length()-2, fname)==0 || l.compare(1, l.length()-2, "GLOBAL")==0) {
				p = 1;
				continue;
			} else {
				p = 0;
			}
		}
		if (p)
		{
			std::size_t sp = l.find('=');
			names.push_back(trim_copy(l.substr(0, sp-1)));
			std::string v = trim_copy(l.substr(sp+1));
			if (v.at(0) == '[' && v.length() == 1) {
				std::vector<std::string> objs;
				while(fin.peek()!=EOF) {
					std::getline(fin, l);
					if (l.at(0) == ']' && l.length() == 1) break;
					objs.push_back(l);
				}

				Array arr = helper::createList(objs.size());
				for(int i = 0; i < objs.size(); i++){
					arr[i] = helper::createString(objs[i].c_str());
				}
				mvalues.push_back(arr);
			} else {
				Value arr; 
				if (std::regex_match(v, std::regex("^[0-9]+(\\.[0-9]+)?$"))) {
					arr = helper::createNumeric(std::stod(v));
				} else {
					arr = helper::createString(v);
				}
				mvalues.push_back(arr);
			}
		}
	}

	NamedArray ret = helper::createNamed(names, mvalues);

	return ret;
}
#if __has_include(<mex.hpp>)
class MexFunction : public matlab::mex::Function {
public:
	std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
  void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
  	matlabPtr = getEngine();

	// Validate arguments
	checkArguments(outputs, inputs);

	// Convert char arrays to string
	matlab::data::StringArray inp = matlabPtr->feval(u"string", 
		  std::vector<matlab::data::Array>({ inputs[0] }));

	//convert from std::basic_string<char16_t> to std::string
	std::wstring_convert<std::codecvt_utf8_utf16<char16_t>,char16_t> convert;

	// Assign outputs
		outputs[0] = parse(convert.to_bytes(inp[0]));
  }

  void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
		if (inputs[0].getType() != matlab::data::ArrayType::CHAR)
		{
		  matlabPtr->feval(u"error", 0, 
			  std::vector<matlab::data::Array>({ factory.createScalar("Input must be a char array") }));
		}

		if (inputs.size() > 1) {
		  matlabPtr->feval(u"error", 0, 
			  std::vector<matlab::data::Array>({ factory.createScalar("parseParams only accepts one input") }));
		}

		if (outputs.size() > 1) {
		  matlabPtr->feval(u"error", 0, 
			  std::vector<matlab::data::Array>({ factory.createScalar("Only one output is returned") }));
		}
  }
};
#elif __has_include("Rcpp.h")
#include "Rcpp.h"

// [[Rcpp::export]]
Rcpp::List parseParams(Rcpp::CharacterVector ff){
	std::string fname = Rcpp::as<std::string>(ff);
	return parse(fname);
}
#else
// bool factory;
#endif