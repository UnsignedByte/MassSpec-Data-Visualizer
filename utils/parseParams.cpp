/*
* @Author: UnsignedByte
* @Date:   03:04:47, 05-Aug-2020
* @Last Modified by:   UnsignedByte
* @Last Modified time: 15:48:58, 26-Aug-2020
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <codecvt>
#include <locale>
#include <regex>

#include <mex.hpp>
#include <mexAdapter.hpp>
// #include <Rcpp.h>

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

// trim from start (copying)
static inline std::string ltrim_copy(std::string s) {
    ltrim(s);
    return s;
}

// trim from end (copying)
static inline std::string rtrim_copy(std::string s) {
    rtrim(s);
    return s;
}

// trim from both ends (copying)
static inline std::string trim_copy(std::string s) {
    trim(s);
    return s;
}

std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
matlab::data::ArrayFactory factory;

struct var {
	int type;
	std::string name;
	double numeric;
	std::string str;
	std::vector<double> numericArr;
	std::vector<std::string> strArr;
};

struct returned {
	matlab::data::Array m;
};

//Hash a cstr
constexpr unsigned int hash(const char *s, int off = 0) {                        
    return !s[off] ? 5381 : (hash(s, off+1)*33) ^ s[off];                           
} 

returned parseParams(const std::string& fname){
	returned ret;

	std::size_t splitter = fname.find_last_of('.');

	std::string type = fname.substr(splitter+1);
	std::string name = fname.substr(0, splitter-1);

	int t;

	switch (hash(type.c_str())) {
		case hash("m"):
			t = 0;
			//ret.m = matlab::data::Struct
			break;
		case hash("r"):
			t = 1;
			break;
		default:
			std::cout << "unsupported type" << std::endl;
			break;
	}

	std::vector<std::string> names;
	std::vector<matlab::data::Array> mvalues;

	std::ifstream fin;
	fin.open("params.p");

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

				switch (t) {
					case 0:
						matlab::data::CellArray arr = factory.createCellArray({1, objs.size()});
						for(int i = 0; i < objs.size(); i++){
							arr[i] = factory.createCharArray(objs[i].c_str());
						}
						mvalues.push_back(arr);
						break;
				}
			} else {
				matlab::data::Array arr; 
				if (std::regex_match(v, std::regex("^[0-9]+(\\.[0-9]+)?$"))) {
					switch (t) {
						case 0:
							arr = factory.createScalar(std::stod(v));
							break;
					}
				} else {
					switch (t) {
						case 0:
							arr = factory.createCharArray(v);
							break;
					}
				}
				mvalues.push_back(arr);
			}
		}
	}

	switch(t) {
		case 0:
			ret.m = factory.createStructArray({1, 1}, names);
			for (int i = 0; i < names.size(); i++){
				ret.m[0][names[i]] = mvalues[i];
			}
			break;
	}

	return ret;
}

class MexFunction : public matlab::mex::Function {
public:
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
    outputs[0] = parseParams(convert.to_bytes(inp[0])).m;
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

