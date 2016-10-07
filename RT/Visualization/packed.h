#ifdef _WIN32 || _WIN64
#define PACKED
#pragma pack(push,1)
#elif __linux__
#define PACKED __attribute__ ((__packed__))
#endif
