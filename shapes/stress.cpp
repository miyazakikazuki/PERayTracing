#include "stress.h"
#include "stats.h"
#include "lowdiscrepancy.h"

namespace pbrt {
class Stress : public Shape {
	Stress(const Transform *ObjectToWorld, const Transform *WorldToObject,
		bool reverseOrientation);
};
}