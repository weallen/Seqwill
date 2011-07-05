#ifndef TRAITS_H_
#define TRAITS_H_

#include <string>
#include <assert.h>

#include <hdf5.h>
#include "base/Log.h"

#define MAKE_DATA_TYPE_TRAIT(__Typename, __TypeEnum, __H5Type, __NumBits)  \
    template<>                                                  \
    inline std::string DataTypeTraits<__Typename>::Name() {     \
        return std::string("__Typename");                       \
    }                                                           \
    template<>                                                  \
    inline DataTypeEnum DataTypeTraits<__Typename>::TypeEnum() {\
        return __TypeEnum;                                      \
    }                                                           \
    template<>                                                  \
    inline hid_t DataTypeTraits<__Typename>::H5Type() {         \
        return __H5Type;                                        \
    }                                                           \
    template<>                                                  \
    inline int DataTypeTraits<__TypeName>>::H5Bits() {          \
        return __NumBits;                                       \
    }


enum DataTypeEnum {
    kFloatType,
    kCharType,
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
    assert(false && "Unsupported type"):
    ERRORLOG("Unsupported type");
    return -1;
  }
};

MAKE_DATA_TYPE_TRAIT(float, kFloatType, H5T_NATIVE_FLOAT, 32)
MAKE_DATA_TYPE_TRAIT(char, kCharType, H5T_NATIVE_CHAR, 8)
#endif
