
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SHAPES_SIMPLEMODEL_H
#define PBRT_SHAPES_SIMPLEMODEL_H

// shapes/simplemodel.h*
#include <map>
#include "shape.h"
#include "stats.h"
#include "triangle.h"

namespace pbrt {

//STAT_MEMORY_COUNTER("Memory/SimpleModel meshes", triMeshBytes);

// Triangle Declarations
class SimpleModel : public Shape {
  public:
    // Triangle Public Methods
    SimpleModel(const Transform *ObjectToWorld, const Transform *WorldToObject,
             bool reverseOrientation, Float width, Float height, Float depth, Float h, 
			int u, int v, int w)
        : Shape(ObjectToWorld, WorldToObject, reverseOrientation), 
        width(width),height(height),depth(depth),
		h(h)
		u(u),v(v), w(w)
    }
    Bounds3f ObjectBound() const;
    Bounds3f WorldBound() const;
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture = true) const;
    bool IntersectP(const Ray &ray, bool testAlphaTexture = true) const;
    Float Area() const;

    using Shape::Sample;  // Bring in the other Sample() overload.
    Interaction Sample(const Point2f &u, Float *pdf) const;

    // Returns the solid angle subtended by the triangle w.r.t. the given
    // reference point p.
    Float SolidAngle(const Point3f &p, int nSamples = 0) const;

	std::vector<std::vector<std::vector<Matrix4x4>>> CalculateStress(), 
  private:
	// Simplemodel Private Data
    const Float width, height, depth;
    const Float h;
    const int u, v, w;
    };

std::vector<std::shared_ptr<Shape>> CreateTriangleMeshShape(
    const Transform *o2w, const Transform *w2o, bool reverseOrientation,
    const ParamSet &params,
    std::map<std::string, std::shared_ptr<Texture<Float>>> *floatTextures =
        nullptr);

}  // namespace pbrt

#endif  // PBRT_SHAPES_SIMPLEMODEL_H
