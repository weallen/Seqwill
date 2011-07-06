#ifndef TRAITS_H_
#define TRAITS_H_

#include <string>
#include <assert.h>

#include <hdf5.h>
#include "base/Log.h"

enum DataTypeEnum {
    kFloatType,
    kCharType,
    kUCharType,
    kPlusMinusStrandType,
    kUnknownType
};

template <typename T>
struct DataTypeTraits {
  static std::string Name() {
    assert(false && "Unsupported type");
    ERRORLOG("Unsupported type");
    return std::string("Unsupported type");
  }

  static DataTypeEnum TypeEnum() {
    assert(false && "Unsupported type");
    ERRORLOG("Unsupported type");
    return kUnknownType;
  }

  static hid_t H5Type() {
    assert(false && "Unsupported type");
    ERRORLOG("Unsupported type");
    return -1;
  }

  static int H5Bits() {
    assert(false && "Unsupported type");
    ERRORLOG("Unsupported type");
    return -1;
  }
};

// Float trait
template<>
inline std::string DataTypeTraits<float>::Name() {
    return std::string("float");
}
template<>
inline DataTypeEnum DataTypeTraits<float>::TypeEnum() {
    return kFloatType;
}
template<>
inline hid_t DataTypeTraits<float>::H5Type() {
    return H5T_NATIVE_FLOAT;
}
template<>
inline int DataTypeTraits<float>::H5Bits() {
    return 32;
}

// Char trait
template<>
inline std::string DataTypeTraits<char>::Name() {
    return std::string("char");
}
template<>
inline DataTypeEnum DataTypeTraits<char>::TypeEnum() {
    return kCharType;
}
template<>
inline hid_t DataTypeTraits<char>::H5Type() {
    return H5T_NATIVE_CHAR;
}
template<>
inline int DataTypeTraits<char>::H5Bits() {
    return 8;
}

// UChar trait
template<>
inline std::string DataTypeTraits<unsigned char>::Name() {
    return std::string("unsigned char");
}
template<>
inline DataTypeEnum DataTypeTraits<unsigned char>::TypeEnum() {
    return kUCharType;
}
template<>
inline hid_t DataTypeTraits<unsigned char>::H5Type() {
    return H5T_NATIVE_UCHAR;
}
template<>
inline int DataTypeTraits<unsigned char>::H5Bits() {
    return 4;
}

#endif
