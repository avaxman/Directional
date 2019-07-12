#ifndef DIRECTIONAL_ASSERT_H
#define DIRECTIONAL_ASSERT_H

#ifndef DIR_ASSERT()
#include <cassert>
#define DIR_ASSERT(x) assert(x)
#endif

#ifndef DIR_ASSERT_M()
#include <cassert>
#define DIR_ASSERT_M(x, m) assert((x) && m)
#endif

#endif