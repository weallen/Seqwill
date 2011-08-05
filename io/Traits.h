#ifndef TRAITS_H_
#define TRAITS_H_

#include <string>
#include <assert.h>

#include <hdf5.h>
#include "base/Log.h"
#include "base/Types.h"

enum DataTypeEnum {
    kFloatType,
    kCharType,
    kUCharType,
    kIntType,
    kPlusMinusDataIntType,
    kPlusMinusDataFloatType,
    kUnknownType
};

template <typename T> struct DataTypeTraits;

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

  static bool IsCompound() {
    assert(false && "Unsupported type");
    ERRORLOG("Unsupported type");
    return false;
  }
};

// plus/minus data trait

template<>
struct DataTypeTraits<PlusMinusDataInt> {
  static std::string Name() {
    return std::string("PlusMinusDataInt");
  }

  static DataTypeEnum TypeEnum() {
    return kPlusMinusDataIntType;
  }

  static hid_t H5Type() {
    hid_t t = H5Tcreate(H5T_COMPOUND, sizeof(PlusMinusDataInt));
    H5Tinsert(t, "plus", HOFFSET(PlusMinusDataInt, plus), H5T_NATIVE_INT);
    H5Tinsert(t, "minus", HOFFSET(PlusMinusDataInt, minus), H5T_NATIVE_INT);
    return t;
  }

  static int H5Bits() {
    return sizeof(PlusMinusDataInt);
  }

  static bool IsCompound() {
    return true;
  }
};

template<>
struct DataTypeTraits<PlusMinusDataFloat> {
  static std::string Name() {
    return std::string("PlusMinusDataFloat");
  }

  static DataTypeEnum TypeEnum() {
    return kPlusMinusDataFloatType;
  }

  static hid_t H5Type() {
    hid_t t = H5Tcreate(H5T_COMPOUND, sizeof(PlusMinusDataFloat));
    H5Tinsert(t, "plus", HOFFSET(PlusMinusDataFloat, plus), H5T_NATIVE_FLOAT);
    H5Tinsert(t, "minus", HOFFSET(PlusMinusDataFloat, minus), H5T_NATIVE_FLOAT);
    return t;
  }

  static int H5Bits() {
    return sizeof(PlusMinusDataFloat);
  }

  static bool IsCompound() {
    return true;
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

template<>
inline bool DataTypeTraits<float>::IsCompound() {
  return false;
}
// int trait
template<>
inline std::string DataTypeTraits<int>::Name() {
  return std::string("int");
}
template<>
inline DataTypeEnum DataTypeTraits<int>::TypeEnum() {
  return kIntType;
}
template<>
inline hid_t DataTypeTraits<int>::H5Type() {
  return H5T_NATIVE_INT;
}
template<>
inline int DataTypeTraits<int>::H5Bits() {
  return sizeof(int);
}

template<>
inline bool DataTypeTraits<int>::IsCompound() {
  return false;
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
template<>
inline bool DataTypeTraits<char>::IsCompound() {
  return false;
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
template<>
inline bool DataTypeTraits<unsigned char>::IsCompound() {
  return false;
}
#endif
