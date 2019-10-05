
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

// shapes/triangle.cpp*
#include "shapes/simplemodel.h"
#include <array>
#include "efloat.h"
#include "paramset.h"
#include "sampling.h"
#include "textures/constant.h"

namespace pbrt {
STAT_PERCENT("Intersections/Ray-triangle intersection tests", nHits, nTests);

// Triangle Method Definitions
STAT_RATIO("Scene/Triangles per triangle mesh", nTris, nMeshes);

Bounds3f SimpleModel::ObjectBound() const {
    // Get triangle vertices in _p0_, _p1_, and _p2_
    return Bounds3f(Point3f(-width / 2.0, -height / 2.0, -depth / 2.0),
                    Point3f(width / 2.0, height / 2.0, depth / 2.0));
}

bool SimpleModel::Intersect(const Ray &ray, Float *tHit,
                            SurfaceInteraction *isect,
                            bool testAlphaTexture) const {
    ProfilePhase p(Prof::TriIntersect);
    ++nTests;
    // Get triangle vertices in _p0_, _p1_, and _p2_
    const Point3f &p0 = mesh->p[_v[0]];
    const Point3f &p1 = mesh->p[_v[1]];
    const Point3f &p2 = mesh->p[_v[2]];

    // Perform ray--triangle intersection test

    // Transform triangle vertices to ray coordinate space

    // Translate vertices based on ray origin
    Point3f p0t = p0 - Vector3f(ray.o);
    Point3f p1t = p1 - Vector3f(ray.o);
    Point3f p2t = p2 - Vector3f(ray.o);

    // Permute components of triangle vertices and ray direction
    int kz = MaxDimension(Abs(ray.d));
    int kx = kz + 1;
    if (kx == 3) kx = 0;
    int ky = kx + 1;
    if (ky == 3) ky = 0;
    Vector3f d = Permute(ray.d, kx, ky, kz);
    p0t = Permute(p0t, kx, ky, kz);
    p1t = Permute(p1t, kx, ky, kz);
    p2t = Permute(p2t, kx, ky, kz);

    // Apply shear transformation to translated vertex positions
    Float Sx = -d.x / d.z;
    Float Sy = -d.y / d.z;
    Float Sz = 1.f / d.z;
    p0t.x += Sx * p0t.z;
    p0t.y += Sy * p0t.z;
    p1t.x += Sx * p1t.z;
    p1t.y += Sy * p1t.z;
    p2t.x += Sx * p2t.z;
    p2t.y += Sy * p2t.z;

    // Compute edge function coefficients _e0_, _e1_, and _e2_
    Float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
    Float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
    Float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

    // Fall back to double precision test at triangle edges
    if (sizeof(Float) == sizeof(float) &&
        (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
        double p2txp1ty = (double)p2t.x * (double)p1t.y;
        double p2typ1tx = (double)p2t.y * (double)p1t.x;
        e0 = (float)(p2typ1tx - p2txp1ty);
        double p0txp2ty = (double)p0t.x * (double)p2t.y;
        double p0typ2tx = (double)p0t.y * (double)p2t.x;
        e1 = (float)(p0typ2tx - p0txp2ty);
        double p1txp0ty = (double)p1t.x * (double)p0t.y;
        double p1typ0tx = (double)p1t.y * (double)p0t.x;
        e2 = (float)(p1typ0tx - p1txp0ty);
    }

    // Perform triangle edge and determinant tests
    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
        return false;
    Float det = e0 + e1 + e2;
    if (det == 0) return false;

    // Compute scaled hit distance to triangle and test against ray $t$ range
    p0t.z *= Sz;
    p1t.z *= Sz;
    p2t.z *= Sz;
    Float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
    if (det < 0 && (tScaled >= 0 || tScaled < ray.tMax * det))
        return false;
    else if (det > 0 && (tScaled <= 0 || tScaled > ray.tMax * det))
        return false;

    // Compute barycentric coordinates and $t$ value for triangle intersection
    Float invDet = 1 / det;
    Float b0 = e0 * invDet;
    Float b1 = e1 * invDet;
    Float b2 = e2 * invDet;
    Float t = tScaled * invDet;

    // Ensure that computed triangle $t$ is conservatively greater than zero

    // Compute $\delta_z$ term for triangle $t$ error bounds
    Float maxZt = MaxComponent(Abs(Vector3f(p0t.z, p1t.z, p2t.z)));
    Float deltaZ = gamma(3) * maxZt;

    // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
    Float maxXt = MaxComponent(Abs(Vector3f(p0t.x, p1t.x, p2t.x)));
    Float maxYt = MaxComponent(Abs(Vector3f(p0t.y, p1t.y, p2t.y)));
    Float deltaX = gamma(5) * (maxXt + maxZt);
    Float deltaY = gamma(5) * (maxYt + maxZt);

    // Compute $\delta_e$ term for triangle $t$ error bounds
    Float deltaE =
        2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

    // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
    Float maxE = MaxComponent(Abs(Vector3f(e0, e1, e2)));
    Float deltaT = 3 *
                   (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
                   std::abs(invDet);
    if (t <= deltaT) return false;

    // Compute triangle partial derivatives
    Vector3f dpdu, dpdv;
    Point2f uv[3];
    GetUVs(uv);

    // Compute deltas for triangle partial derivatives
    Vector2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
    Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
    Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
    bool degenerateUV = std::abs(determinant) < 1e-8;
    if (!degenerateUV) {
        Float invdet = 1 / determinant;
        dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
        dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
    }
    if (degenerateUV || Cross(dpdu, dpdv).LengthSquared() == 0) {
        // Handle zero determinant for triangle partial derivative matrix
        Vector3f ng = Cross(p2 - p0, p1 - p0);
        if (ng.LengthSquared() == 0)
            // The triangle is actually degenerate; the intersection is
            // bogus.
            return false;

        CoordinateSystem(Normalize(ng), &dpdu, &dpdv);
    }

    // Compute error bounds for triangle intersection
    Float xAbsSum =
        (std::abs(b0 * p0.x) + std::abs(b1 * p1.x) + std::abs(b2 * p2.x));
    Float yAbsSum =
        (std::abs(b0 * p0.y) + std::abs(b1 * p1.y) + std::abs(b2 * p2.y));
    Float zAbsSum =
        (std::abs(b0 * p0.z) + std::abs(b1 * p1.z) + std::abs(b2 * p2.z));
    Vector3f pError = gamma(7) * Vector3f(xAbsSum, yAbsSum, zAbsSum);

    // Interpolate $(u,v)$ parametric coordinates and hit point
    Point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
    Point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

    // Test intersection against alpha texture, if present
    if (testAlphaTexture && mesh->alphaMask) {
        SurfaceInteraction isectLocal(pHit, Vector3f(0, 0, 0), uvHit, -ray.d,
                                      dpdu, dpdv, Normal3f(0, 0, 0),
                                      Normal3f(0, 0, 0), ray.time, this);
        if (mesh->alphaMask->Evaluate(isectLocal) == 0) return false;
    }

    // Fill in _SurfaceInteraction_ from triangle hit
    *isect = SurfaceInteraction(pHit, pError, uvHit, -ray.d, dpdu, dpdv,
                                Normal3f(0, 0, 0), Normal3f(0, 0, 0), ray.time,
                                this, faceIndex);

    // Override surface normal in _isect_ for triangle
    isect->n = isect->shading.n = Normal3f(Normalize(Cross(dp02, dp12)));
    if (mesh->n || mesh->s) {
        // Initialize _Triangle_ shading geometry

        // Compute shading normal _ns_ for triangle
        Normal3f ns;
        if (mesh->n) {
            ns = (b0 * mesh->n[_v[0]] + b1 * mesh->n[_v[1]] +
                  b2 * mesh->n[_v[2]]);
            if (ns.LengthSquared() > 0)
                ns = Normalize(ns);
            else
                ns = isect->n;
        } else
            ns = isect->n;

        // Compute shading tangent _ss_ for triangle
        Vector3f ss;
        if (mesh->s) {
            ss = (b0 * mesh->s[_v[0]] + b1 * mesh->s[_v[1]] +
                  b2 * mesh->s[_v[2]]);
            if (ss.LengthSquared() > 0)
                ss = Normalize(ss);
            else
                ss = Normalize(isect->dpdu);
        } else
            ss = Normalize(isect->dpdu);

        // Compute shading bitangent _ts_ for triangle and adjust _ss_
        Vector3f ts = Cross(ss, ns);
        if (ts.LengthSquared() > 0.f) {
            ts = Normalize(ts);
            ss = Cross(ts, ns);
        } else
            CoordinateSystem((Vector3f)ns, &ss, &ts);

        // Compute $\dndu$ and $\dndv$ for triangle shading geometry
        Normal3f dndu, dndv;
        if (mesh->n) {
            // Compute deltas for triangle partial derivatives of normal
            Vector2f duv02 = uv[0] - uv[2];
            Vector2f duv12 = uv[1] - uv[2];
            Normal3f dn1 = mesh->n[_v[0]] - mesh->n[_v[2]];
            Normal3f dn2 = mesh->n[_v[1]] - mesh->n[_v[2]];
            Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
            bool degenerateUV = std::abs(determinant) < 1e-8;
            if (degenerateUV) {
                // We can still compute dndu and dndv, with respect to the
                // same arbitrary coordinate system we use to compute dpdu
                // and dpdv when this happens. It's important to do this
                // (rather than giving up) so that ray differentials for
                // rays reflected from triangles with degenerate
                // parameterizations are still reasonable.
                Vector3f dn = Cross(Vector3f(mesh->n[_v[2]] - mesh->n[_v[0]]),
                                    Vector3f(mesh->n[_v[1]] - mesh->n[_v[0]]));
                if (dn.LengthSquared() == 0)
                    dndu = dndv = Normal3f(0, 0, 0);
                else {
                    Vector3f dnu, dnv;
                    CoordinateSystem(dn, &dnu, &dnv);
                    dndu = Normal3f(dnu);
                    dndv = Normal3f(dnv);
                }
            } else {
                Float invDet = 1 / determinant;
                dndu = (duv12[1] * dn1 - duv02[1] * dn2) * invDet;
                dndv = (-duv12[0] * dn1 + duv02[0] * dn2) * invDet;
            }
        } else
            dndu = dndv = Normal3f(0, 0, 0);
        isect->SetShadingGeometry(ss, ts, dndu, dndv, true);
    }

    // Ensure correct orientation of the geometric normal
    if (mesh->n)
        isect->n = Faceforward(isect->n, isect->shading.n);
    else if (reverseOrientation ^ transformSwapsHandedness)
        isect->n = isect->shading.n = -isect->n;
    *tHit = t;
    ++nHits;
    return true;
}

bool SimpleModel::IntersectP(const Ray &ray, bool testAlphaTexture) const {
    ProfilePhase p(Prof::TriIntersectP);
    ++nTests;
    // Get triangle vertices in _p0_, _p1_, and _p2_
    const Point3f &p0 = mesh->p[_v[0]];
    const Point3f &p1 = mesh->p[_v[1]];
    const Point3f &p2 = mesh->p[_v[2]];

    // Perform ray--triangle intersection test

    // Transform triangle vertices to ray coordinate space

    // Translate vertices based on ray origin
    Point3f p0t = p0 - Vector3f(ray.o);
    Point3f p1t = p1 - Vector3f(ray.o);
    Point3f p2t = p2 - Vector3f(ray.o);

    // Permute components of triangle vertices and ray direction
    int kz = MaxDimension(Abs(ray.d));
    int kx = kz + 1;
    if (kx == 3) kx = 0;
    int ky = kx + 1;
    if (ky == 3) ky = 0;
    Vector3f d = Permute(ray.d, kx, ky, kz);
    p0t = Permute(p0t, kx, ky, kz);
    p1t = Permute(p1t, kx, ky, kz);
    p2t = Permute(p2t, kx, ky, kz);

    // Apply shear transformation to translated vertex positions
    Float Sx = -d.x / d.z;
    Float Sy = -d.y / d.z;
    Float Sz = 1.f / d.z;
    p0t.x += Sx * p0t.z;
    p0t.y += Sy * p0t.z;
    p1t.x += Sx * p1t.z;
    p1t.y += Sy * p1t.z;
    p2t.x += Sx * p2t.z;
    p2t.y += Sy * p2t.z;

    // Compute edge function coefficients _e0_, _e1_, and _e2_
    Float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
    Float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
    Float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

    // Fall back to double precision test at triangle edges
    if (sizeof(Float) == sizeof(float) &&
        (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
        double p2txp1ty = (double)p2t.x * (double)p1t.y;
        double p2typ1tx = (double)p2t.y * (double)p1t.x;
        e0 = (float)(p2typ1tx - p2txp1ty);
        double p0txp2ty = (double)p0t.x * (double)p2t.y;
        double p0typ2tx = (double)p0t.y * (double)p2t.x;
        e1 = (float)(p0typ2tx - p0txp2ty);
        double p1txp0ty = (double)p1t.x * (double)p0t.y;
        double p1typ0tx = (double)p1t.y * (double)p0t.x;
        e2 = (float)(p1typ0tx - p1txp0ty);
    }

    // Perform triangle edge and determinant tests
    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
        return false;
    Float det = e0 + e1 + e2;
    if (det == 0) return false;

    // Compute scaled hit distance to triangle and test against ray $t$ range
    p0t.z *= Sz;
    p1t.z *= Sz;
    p2t.z *= Sz;
    Float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
    if (det < 0 && (tScaled >= 0 || tScaled < ray.tMax * det))
        return false;
    else if (det > 0 && (tScaled <= 0 || tScaled > ray.tMax * det))
        return false;

    // Compute barycentric coordinates and $t$ value for triangle intersection
    Float invDet = 1 / det;
    Float b0 = e0 * invDet;
    Float b1 = e1 * invDet;
    Float b2 = e2 * invDet;
    Float t = tScaled * invDet;

    // Ensure that computed triangle $t$ is conservatively greater than zero

    // Compute $\delta_z$ term for triangle $t$ error bounds
    Float maxZt = MaxComponent(Abs(Vector3f(p0t.z, p1t.z, p2t.z)));
    Float deltaZ = gamma(3) * maxZt;

    // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
    Float maxXt = MaxComponent(Abs(Vector3f(p0t.x, p1t.x, p2t.x)));
    Float maxYt = MaxComponent(Abs(Vector3f(p0t.y, p1t.y, p2t.y)));
    Float deltaX = gamma(5) * (maxXt + maxZt);
    Float deltaY = gamma(5) * (maxYt + maxZt);

    // Compute $\delta_e$ term for triangle $t$ error bounds
    Float deltaE =
        2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

    // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
    Float maxE = MaxComponent(Abs(Vector3f(e0, e1, e2)));
    Float deltaT = 3 *
                   (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
                   std::abs(invDet);
    if (t <= deltaT) return false;

    // Test shadow ray intersection against alpha texture, if present
    if (testAlphaTexture && (mesh->alphaMask || mesh->shadowAlphaMask)) {
        // Compute triangle partial derivatives
        Vector3f dpdu, dpdv;
        Point2f uv[3];
        GetUVs(uv);

        // Compute deltas for triangle partial derivatives
        Vector2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
        Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
        Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
        bool degenerateUV = std::abs(determinant) < 1e-8;
        if (!degenerateUV) {
            Float invdet = 1 / determinant;
            dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
            dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
        }
        if (degenerateUV || Cross(dpdu, dpdv).LengthSquared() == 0) {
            // Handle zero determinant for triangle partial derivative matrix
            Vector3f ng = Cross(p2 - p0, p1 - p0);
            if (ng.LengthSquared() == 0)
                // The triangle is actually degenerate; the intersection is
                // bogus.
                return false;

            CoordinateSystem(Normalize(Cross(p2 - p0, p1 - p0)), &dpdu, &dpdv);
        }

        // Interpolate $(u,v)$ parametric coordinates and hit point
        Point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
        Point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];
        SurfaceInteraction isectLocal(pHit, Vector3f(0, 0, 0), uvHit, -ray.d,
                                      dpdu, dpdv, Normal3f(0, 0, 0),
                                      Normal3f(0, 0, 0), ray.time, this);
        if (mesh->alphaMask && mesh->alphaMask->Evaluate(isectLocal) == 0)
            return false;
        if (mesh->shadowAlphaMask &&
            mesh->shadowAlphaMask->Evaluate(isectLocal) == 0)
            return false;
    }
    ++nHits;
    return true;
}

Float SimpleModel::Area() const {
    const Point3f &p0 = mesh->p[_v[0]];
    const Point3f &p1 = mesh->p[_v[1]];
    const Point3f &p2 = mesh->p[_v[2]];
    return 0.5 * Cross(p1 - p0, p2 - p0).Length();
}

Interaction SimpleModel::Sample(const Point2f &u, Float *pdf) const {
    Point2f b = UniformSampleTriangle(u);
    // Get triangle vertices in _p0_, _p1_, and _p2_
    const Point3f &p0 = mesh->p[_v[0]];
    const Point3f &p1 = mesh->p[_v[1]];
    const Point3f &p2 = mesh->p[_v[2]];
    Interaction it;
    it.p = b[0] * p0 + b[1] * p1 + (1 - b[0] - b[1]) * p2;
    // Compute surface normal for sampled point on triangle
    it.n = Normalize(Normal3f(Cross(p1 - p0, p2 - p0)));
    // Ensure correct orientation of the geometric normal; follow the same
    // approach as was used in SimpleModel::Intersect().
    if (mesh->n) {
        Normal3f ns(b[0] * mesh->n[_v[0]] + b[1] * mesh->n[_v[1]] +
                    (1 - b[0] - b[1]) * mesh->n[_v[2]]);
        it.n = Faceforward(it.n, ns);
    } else if (reverseOrientation ^ transformSwapsHandedness)
        it.n *= -1;

    // Compute error bounds for sampled point on triangle
    Point3f pAbsSum =
        Abs(b[0] * p0) + Abs(b[1] * p1) + Abs((1 - b[0] - b[1]) * p2);
    it.pError = gamma(6) * Vector3f(pAbsSum.x, pAbsSum.y, pAbsSum.z);
    *pdf = 1 / Area();
    return it;
}
/*
Float SimpleModel::SolidAngle(const Point3f &p, int nSamples) const {
    // Project the vertices into the unit sphere around p.
    std::array<Vector3f, 3> pSphere = {Normalize(mesh->p[v[0]] - p),
                                       Normalize(mesh->p[v[1]] - p),
                                       Normalize(mesh->p[v[2]] - p)};

    // http://math.stackexchange.com/questions/9819/area-of-a-spherical-triangle
    // Girard's theorem: surface area of a spherical triangle on a unit
    // sphere is the 'excess angle' alpha+beta+gamma-pi, where
    // alpha/beta/gamma are the interior angles at the vertices.
    //
    // Given three vertices on the sphere, a, b, c, then we can compute,
    // for example, the angle c->a->b by
    //
    // cos theta =  Dot(Cross(c, a), Cross(b, a)) /
    //              (Length(Cross(c, a)) * Length(Cross(b, a))).
    //
    Vector3f cross01 = (Cross(pSphere[0], pSphere[1]));
    Vector3f cross12 = (Cross(pSphere[1], pSphere[2]));
    Vector3f cross20 = (Cross(pSphere[2], pSphere[0]));

    // Some of these vectors may be degenerate. In this case, we don't want
    // to normalize them so that we don't hit an assert. This is fine,
    // since the corresponding dot products below will be zero.
    if (cross01.LengthSquared() > 0) cross01 = Normalize(cross01);
    if (cross12.LengthSquared() > 0) cross12 = Normalize(cross12);
    if (cross20.LengthSquared() > 0) cross20 = Normalize(cross20);

    // We only need to do three cross products to evaluate the angles at
    // all three vertices, though, since we can take advantage of the fact
    // that Cross(a, b) = -Cross(b, a).
    return std::abs(std::acos(Clamp(Dot(cross01, -cross12), -1, 1)) +
                    std::acos(Clamp(Dot(cross12, -cross20), -1, 1)) +
                    std::acos(Clamp(Dot(cross20, -cross01), -1, 1)) - Pi);
}*/
std::vector<std::vector<std::vector<Matrix4x4>>> SimpleModel::CalculateStress(
    const int u, const int v, const int w) {
    std::vector<Point3f> vertices;
    for (int i = 0; i < u; i++) {
        for (int j = 0; j < v; j++) {
            for (int k = 0; k < w; k++) {
                vertices.push_back(Point3f(i, j, k));
            }
        }
    }

    /*foreach(Vector3 vertex in vertices) {
            Debug.Log(vertex.x);
            Debug.Log(vertex.y);
            Debug.Log(vertex.z);
    }*/

    std::vector<std::vector<int>> indices;
    std::vector<int> verticesIndices;
    int kuv, uv = u * v, ju;
    for (int i = 0; i < u; i++) {
        for (int j = 0; j < v; j++) {
            ju = j * u;
            for (int k = 0; k < w; k++, verticesIndices = {}) {
                kuv = k * u * v;
                verticesIndices.push_back(kuv + ju + i);
                verticesIndices.push_back(kuv + ju + i + 1);
                verticesIndices.push_back(kuv + ju + i + u);
                verticesIndices.push_back(kuv + ju + i + u + 1);
                verticesIndices.push_back(kuv + ju + i + uv);
                verticesIndices.push_back(kuv + ju + i + uv + 1);
                verticesIndices.push_back(kuv + ju + i + uv + u);
                verticesIndices.push_back(kuv + ju + i + uv + u + 1);
            }
            indices.push_back(verticesIndices);
        }
    }

    /*
    foreach (int[] index in indices) {
            Debug.Log(index[0]);
            Debug.Log(index[1]);
            Debug.Log(index[2]);
    }*/
    Float a = -sqrtf(1.0f / 3.0f), b = sqrtf(1.0f / 3.0f);

    std::vector<Vector3f> r = {Vector3f(a, a, a), Vector3f(a, b, a),
                               Vector3f(b, a, a), Vector3f(b, b, a),
                               Vector3f(a, a, b), Vector3f(a, b, b),
                               Vector3f(b, a, b), Vector3f(b, b, b)};
    int E = 100;     // ÉÑÉìÉOó¶
    Float Nu = 0.5;  // É|ÉAÉ\Éìî‰
    Float _Nu = 1 - Nu * Nu;
    Float G = E / (2 * (1.0 - Nu));  // ÇπÇÒífíeê´åWêî
    Float delta = 0.5;

    Float onemNu = 1 - Nu;
    Float onem2Nu = 1 - 2 * Nu;
    Float onepNu = 1 + Nu;
    Float EdopNuom2Nu = E / (onepNu * onem2Nu);
    Float c = onemNu * EdopNuom2Nu;
    Float d = Nu * EdopNuom2Nu;
    Float e = onem2Nu / 2.0 * EdopNuom2Nu;

    std::vector<std::vector<Float>> D = {
        {c, d, d, 0.0, 0.0, 0.0},     {d, c, d, 0.0, 0.0, 0.0},
        {d, d, c, 0.0, 0.0, 0.0},     {0.0, 0.0, 0.0, e, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, e, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, e}};

    std::vector<Float> F(3 * u * v * w);
    std::vector<std::vector<std::vector<Matrix4x4>>> result;
    std::vector<Float> N, dNdr1, dNdr2, dNdr3;
    Float r1, r2, r3, r1p, r2p, r3p, r1m, r2m, r3m;
    Float xg[2] = {-1.0f / sqrtf(3.0), 1.0 / sqrtf(3.0)};
    for (std::vector<int> index : indices) {
        for (int i = 0; i < 2; i++) {
            r1 = xg[i];
            for (int j = 0; j < 2; i++) {
                r2 = xg[j];
                for (int _k = 0; _k < 2; _k++) {
                    r3 = xg[_k];
                    r1p = 1 + r1;
                    r2p = 1 + r2;
                    r3p = 1 + r3;
                    r1m = 1 - r1;
                    r2m = 1 - r2;
                    r3m = 1 - r3;
                    // interpolation function
                    N[0] = r1m * r2m * r3m / 8.0;
                    N[1] = r1p * r2m * r3m / 8.0;
                    N[2] = r1p * r2p * r3m / 8.0;
                    N[3] = r1m * r2p * r3m / 8.0;
                    N[4] = r1m * r2m * r3p / 8.0;
                    N[5] = r1p * r2m * r3p / 8.0;
                    N[6] = r1p * r2p * r3p / 8.0;
                    N[7] = r1m * r2p * r3p / 8.0;
                    // derivative of interpolation function for r1
                    dNdr1[0] = -r2m * r3m / 8.0;
                    dNdr1[1] = r2m * r3m / 8.0;
                    dNdr1[2] = r2p * r3m / 8.0;
                    dNdr1[3] = -r2p * r3m / 8.0;
                    dNdr1[4] = -r2m * r3p / 8.0;
                    dNdr1[5] = r2m * r3p / 8.0;
                    dNdr1[6] = r2p * r3p / 8.0;
                    dNdr1[7] = -r2p * r3p / 8.0;
                    // derivative of interpolation function for r2
                    dNdr2[0] = -r1m * r3m / 8.0;
                    dNdr2[1] = -r1m * r3m / 8.0;
                    dNdr2[2] = r1p * r3m / 8.0;
                    dNdr2[3] = r1p * r3m / 8.0;
                    dNdr2[4] = -r1m * r3p / 8.0;
                    dNdr2[5] = -r1m * r3p / 8.0;
                    dNdr2[6] = r1p * r3p / 8.0;
                    dNdr2[7] = r1p * r3p / 8.0;
                    // derivative of interpolation function for r3
                    dNdr1[0] = -r2m * r1m / 8.0;
                    dNdr1[1] = -r2m * r1m / 8.0;
                    dNdr1[2] = -r2p * r1m / 8.0;
                    dNdr1[3] = -r2p * r1m / 8.0;
                    dNdr1[4] = r2m * r1p / 8.0;
                    dNdr1[5] = r2m * r1p / 8.0;
                    dNdr1[6] = r2p * r1p / 8.0;
                    dNdr1[7] = r2p * r1p / 8.0;
                    // Jacobi Matrix
                    std::vector<std::vector<Float>> J = {
                        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
                    for (int l = 0; l < 8; l++) {
                        J[0][0] += dNdr1[l] * vertices[index[l]][0];
                        J[0][1] += dNdr1[l] * vertices[index[l]][1];
                        J[0][2] += dNdr1[l] * vertices[index[l]][2];
                        J[1][0] += dNdr1[l] * vertices[index[l]][0];
                        J[1][1] += dNdr1[l] * vertices[index[l]][1];
                        J[1][2] += dNdr1[l] * vertices[index[l]][2];
                        J[2][0] += dNdr1[l] * vertices[index[l]][0];
                        J[2][1] += dNdr1[l] * vertices[index[l]][1];
                        J[2][2] += dNdr1[l] * vertices[index[l]][2];
                    }
                    // inverse matrix to J
                    std::vector<std::vector<Float>> invJ;

                    for (int l = 0; l < 3; l++) {
                        for (int m = 0; m < 3; m++) {
                            invJ[l][m] = (l == m) ? 1.0 : 0.0;
                        }
                    }
                    Float buf;
                    std::vector<std::vector<Float>> tmpJ = {
                        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
                    for (int l = 0; l < 3; l++) {
                        for (int m = 0; m < 3; m++) {
                            tmpJ[l][m] = J[l][m];
                        }
                    }

                    for (int l = 0; l < 3; l++) {
                        buf = 1 / tmpJ[l][l];
                        for (int m = 0; m < 3; m++) {
                            tmpJ[l][m] *= buf;
                            invJ[l][m] *= buf;
                        }
                        for (int m = 0; m < 3; m++) {
                            if (l != m) {
                                buf = tmpJ[m][l];
                                for (int n = 0; n < 3; n++) {
                                    tmpJ[m][n] -= tmpJ[l][n] * buf;
                                    invJ[m][n] -= invJ[l][n] * buf;
                                }
                            }
                        }
                    }
                    // determinant of J
                    Float det = 0.0;
                    det = J[0][0] * J[1][1] * J[2][2] +
                          J[1][0] * J[2][1] * J[0][2] +
                          J[2][0] * J[0][1] * J[1][2] -
                          J[2][0] * J[1][1] * J[0][2] -
                          J[1][0] * J[0][1] * J[2][2] -
                          J[0][0] * J[2][1] * J[1][2];

                    // B and its transpose matrix
                    std::vector<std::vector<Float>> B(6,
                                                      std::vector<Float>(24));
                    std::vector<std::vector<Float>> BT(24,
                                                       std::vector<Float>(6));

                    for (int l = 0; l < 8; l++) {
                        Float dNdx = invJ[0][0] * dNdr1[l] +
                                     invJ[0][1] * dNdr2[l] + invJ[0][2];
                        Float dNdy = invJ[1][0] * dNdr1[l] +
                                     invJ[1][1] * dNdr2[l] + invJ[1][2];
                        Float dNdz = invJ[2][0] * dNdr1[l] +
                                     invJ[2][1] * dNdr2[l] + invJ[2][2];

                        B[0][3 * l] = dNdx;
                        B[1][3 * l] = 0.0;
                        B[2][3 * l] = 0.0;
                        B[3][3 * l] = dNdy;
                        B[4][3 * l] = 0.0;
                        B[5][3 * l] = dNdz;
                        B[0][3 * l + 1] = 0.0;
                        B[1][3 * l + 1] = dNdy;
                        B[2][3 * l + 1] = 0.0;
                        B[3][3 * l + 1] = dNdx;
                        B[4][3 * l + 1] = dNdz;
                        B[5][3 * l + 1] = 0.0;
                        B[0][3 * l + 2] = 0.0;
                        B[1][3 * l + 2] = 0.0;
                        B[2][3 * l + 2] = dNdz;
                        B[3][3 * l + 2] = 0.0;
                        B[4][3 * l + 2] = dNdy;
                        B[5][3 * l + 2] = dNdx;

                        BT[3 * l][0] = dNdx;
                        BT[3 * l][1] = 0.0;
                        BT[3 * l][2] = 0.0;
                        BT[3 * l][3] = dNdy;
                        BT[3 * l][4] = 0.0;
                        BT[3 * l][5] = dNdz;
                        BT[3 * l + 1][0] = 0.0;
                        BT[3 * l + 1][1] = dNdy;
                        BT[3 * l + 1][2] = 0.0;
                        BT[3 * l + 1][3] = dNdx;
                        BT[3 * l + 1][4] = dNdz;
                        BT[3 * l + 1][5] = 0.0;
                        BT[3 * l + 2][0] = 0.0;
                        BT[3 * l + 2][1] = 0.0;
                        BT[3 * l + 2][2] = dNdz;
                        BT[3 * l + 2][3] = 0.0;
                        BT[3 * l + 2][4] = dNdy;
                        BT[3 * l + 2][5] = dNdx;
                    }
                    // stiffness matrix
                    std::vector<std::vector<Float>> tmp(6,
                                                        std::vector<Float>(6));
                    std::vector<std::vector<Float>> k(24,
                                                      std::vector<Float>(24));
                    std::vector<std::vector<Float>> K(
                        3 * u * v * w, std::vector<Float>(3 * u * v * w));

                    for (int l = 0; l < 24; l++) {
                        for (int m = 0; m < 6; m++) {
                            for (int n = 0; n < 6; n++) {
                                tmp[l][m] += det * BT[l][n] * D[n][m];
                            }
                        }
                    }
                    for (int l = 0; l < 24; l++) {
                        for (int m = 0; m < 24; m++) {
                            for (int n = 0; n < 6; n++) {
                                k[l][m] += tmp[l][n] * B[n][m];
                            }
                        }
                    }

                    for (int l = 0; l < 24; l++) {
                        for (int m = 0; m < 24; m++) {
                            K[8 * index[l / 3] + l % 3]
                             [8 * index[m / 3] + m % 3] = k[l][m];
                        }
                    }
                }
            }
        }
    }

    // ltype = 1  BX: BODY FORCE IN X-DIRECTION
    // ltype = 2  BY: BODY FORCE IN Y-DIRECTION
    // ltype = 3  BZ: BODY FORCE IN Z-DIRECTION
    // ltype = 10 P1: TRACTIOM IN NORMAL-DIRECTION FOR FACE-1
    // ltype = 20 P2: TRACTIOM IN NORMAL-DIRECTION FOR FACE-2
    // ltype = 30 P3: TRACTIOM IN NORMAL-DIRECTION FOR FACE-3
    // ltype = 40 P4: TRACTIOM IN NORMAL-DIRECTION FOR FACE-4
    // ltype = 50 P5: TRACTIOM IN NORMAL-DIRECTION FOR FACE-5
    // ltype = 60 P6: TRACTIOM IN NORMAL-DIRECTION FOR FACE-6
    for (std::vector<int> index : indices) {
        int ltype = 10;
        bool issuf = false;
        std::vector<int> nL(4), nI(4);
        if (ltype < 10) {
            issuf = false;
        } else {
            issuf = true;
        };

        if (issuf) {
            switch (ltype) {
            case 10:
                nL = {index[0], index[1], index[2], index[3]};
                nI = {0, 1, 2, 3};
            case 20:
                nL = {index[7], index[6], index[5], index[4]};
                nI = {7, 6, 5, 4};
            case 30:
                nL = {index[4], index[5], index[1], index[0]};
                nI = {4, 5, 1, 0};
            case 40:
                nL = {index[5], index[6], index[2], index[1]};
                nI = {5, 6, 2, 1};
            case 50:
                nL = {index[6], index[7], index[3], index[2]};
                nI = {6, 7, 3, 2};
            case 60:
                nL = {index[7], index[4], index[0], index[3]};
                nI = {7, 4, 0, 3};
            default:
                break;
            }
            for (int i = 0; i < 2; i++) {
                r1 = xg[i];
                for (int j = 0; j < 2; i++) {
                    r2 = xg[j];
                    r1p = 1 + r1;
                    r2p = 1 + r2;
                    r1m = 1 - r1;
                    r2m = 1 - r2;
                    // interpolation function
                    N[0] = r1m * r2m / 4.0;
                    N[1] = r1p * r2m / 4.0;
                    N[2] = r1p * r2p / 4.0;
                    N[3] = r1m * r2p / 4.0;
                    // derivative of interpolation function for r1
                    dNdr1[0] = -r2m / 4.0;
                    dNdr1[1] = r2m / 4.0;
                    dNdr1[2] = r2p / 4.0;
                    dNdr1[3] = -r2p / 4.0;
                    // derivative of interpolation function for r2
                    dNdr2[0] = -r1m / 4.0;
                    dNdr2[1] = -r1m / 4.0;
                    dNdr2[2] = r1p / 4.0;
                    dNdr2[3] = r1p / 4.0;
                    // G1,G2
                    std::vector<Float> g1 = {0.0, 0.0, 0.0};
                    std::vector<Float> g2 = {0.0, 0.0, 0.0};
                    for (int _k = 0; _k < 4; _k++) {
                        for (int l = 0; l < 3; l++) {
                            g1[l] += dNdr1[_k] * vertices[nL[_k]][l];
                            g2[l] += dNdr2[_k] * vertices[nL[_k]][l];
                        }
                    }
                    std::vector<Float> g3 = {g1[1] * g2[2] - g1[2] * g2[1],
                                             g1[2] * g2[0] - g1[0] * g2[2],
                                             g1[0] * g2[1] - g1[1] * g2[0]};
                    // Jacobi Matrix
                    std::vector<std::vector<Float>> J = {g1, g2, g3};
                    Float det = 0.0;
                    det = J[0][0] * J[1][1] * J[2][2] +
                          J[1][0] * J[2][1] * J[0][2] +
                          J[2][0] * J[0][1] * J[1][2] -
                          J[2][0] * J[1][1] * J[0][2] -
                          J[1][0] * J[0][1] * J[2][2] -
                          J[0][0] * J[2][1] * J[1][2];

                    Float weight = det;
                    for (int l = 0; l < 4; l++) {
                        F[3 * nL[l]] += val /*strength of traction*/ * weight *
                                        N[l] * g3[0];
                        F[3 * nL[l] + 1] += val /*strength of traction*/ *
                                            weight * N[l] * g3[1];
                        F[3 * nL[l] + 2] += val /*strength of traction*/ *
                                            weight * N[l] * g3[2];
                    }
                }
            }
        }
    }
    // Ignore Body Force

    // Gaussian Elimination
    int n = 3 * u * v * w;
    for (int k = 0; k < n; ++k) {
        double Kkk = K[k][k];
        for (int i = k + 1; i < n; ++i) {
            double Kik = K[i][k];
            for (int j = k; j < n + 1; ++j) {
                K[i][j] = K[i][j] - Kik * (K[k][j] / Kkk);
            }
            F[i] = F[i] - Kik * (F[i] / Kkk);
        }
    }
    F[n] = F[n] / K[n - 1][n - 1];
    for (int i = n - 2; i >= 0; --i) {
        double ax = 0.0;
        for (int j = i + 1; j < n; ++j) {
            ax += K[i][j] * F[j];
        }
        F[i] = (F[i] - ax) / K[i][i];
    }

    for (std::vector<int> index : indices) {
        for (int i = 0; i < 2; i++) {
            r1 = xg[i];
            for (int j = 0; j < 2; i++) {
                r2 = xg[j];
                for (int _k = 0; _k < 2; _k++) {
                    r3 = xg[_k];
                    r1p = 1 + r1;
                    r2p = 1 + r2;
                    r3p = 1 + r3;
                    r1m = 1 - r1;
                    r2m = 1 - r2;
                    r3m = 1 - r3;
                    // interpolation function
                    N[0] = r1m * r2m * r3m / 8.0;
                    N[1] = r1p * r2m * r3m / 8.0;
                    N[2] = r1p * r2p * r3m / 8.0;
                    N[3] = r1m * r2p * r3m / 8.0;
                    N[4] = r1m * r2m * r3p / 8.0;
                    N[5] = r1p * r2m * r3p / 8.0;
                    N[6] = r1p * r2p * r3p / 8.0;
                    N[7] = r1m * r2p * r3p / 8.0;
                    // derivative of interpolation function for r1
                    dNdr1[0] = -r2m * r3m / 8.0;
                    dNdr1[1] = r2m * r3m / 8.0;
                    dNdr1[2] = r2p * r3m / 8.0;
                    dNdr1[3] = -r2p * r3m / 8.0;
                    dNdr1[4] = -r2m * r3p / 8.0;
                    dNdr1[5] = r2m * r3p / 8.0;
                    dNdr1[6] = r2p * r3p / 8.0;
                    dNdr1[7] = -r2p * r3p / 8.0;
                    // derivative of interpolation function for r2
                    dNdr2[0] = -r1m * r3m / 8.0;
                    dNdr2[1] = -r1m * r3m / 8.0;
                    dNdr2[2] = r1p * r3m / 8.0;
                    dNdr2[3] = r1p * r3m / 8.0;
                    dNdr2[4] = -r1m * r3p / 8.0;
                    dNdr2[5] = -r1m * r3p / 8.0;
                    dNdr2[6] = r1p * r3p / 8.0;
                    dNdr2[7] = r1p * r3p / 8.0;
                    // derivative of interpolation function for r3
                    dNdr1[0] = -r2m * r1m / 8.0;
                    dNdr1[1] = -r2m * r1m / 8.0;
                    dNdr1[2] = -r2p * r1m / 8.0;
                    dNdr1[3] = -r2p * r1m / 8.0;
                    dNdr1[4] = r2m * r1p / 8.0;
                    dNdr1[5] = r2m * r1p / 8.0;
                    dNdr1[6] = r2p * r1p / 8.0;
                    dNdr1[7] = r2p * r1p / 8.0;
                    // Jacobi Matrix
                    std::vector<std::vector<Float>> J = {
                        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
                    for (int l = 0; l < 8; l++) {
                        J[0][0] += dNdr1[l] * vertices[index[l]][0];
                        J[0][1] += dNdr1[l] * vertices[index[l]][1];
                        J[0][2] += dNdr1[l] * vertices[index[l]][2];
                        J[1][0] += dNdr1[l] * vertices[index[l]][0];
                        J[1][1] += dNdr1[l] * vertices[index[l]][1];
                        J[1][2] += dNdr1[l] * vertices[index[l]][2];
                        J[2][0] += dNdr1[l] * vertices[index[l]][0];
                        J[2][1] += dNdr1[l] * vertices[index[l]][1];
                        J[2][2] += dNdr1[l] * vertices[index[l]][2];
                    }
                    // inverse matrix to J
                    std::vector<std::vector<Float>> invJ;

                    for (int l = 0; l < 3; l++) {
                        for (int m = 0; m < 3; m++) {
                            invJ[l][m] = (l == m) ? 1.0 : 0.0;
                        }
                    }
                    Float buf;
                    std::vector<std::vector<Float>> tmpJ = {
                        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
                    for (int l = 0; l < 3; l++) {
                        for (int m = 0; m < 3; m++) {
                            tmpJ[l][m] = J[l][m];
                        }
                    }

                    for (int l = 0; l < 3; l++) {
                        buf = 1 / tmpJ[l][l];
                        for (int m = 0; m < 3; m++) {
                            tmpJ[l][m] *= buf;
                            invJ[l][m] *= buf;
                        }
                        for (int m = 0; m < 3; m++) {
                            if (l != m) {
                                buf = tmpJ[m][l];
                                for (int n = 0; n < 3; n++) {
                                    tmpJ[m][n] -= tmpJ[l][n] * buf;
                                    invJ[m][n] -= invJ[l][n] * buf;
                                }
                            }
                        }
                    }
                    // determinant of J
                    Float det = 0.0;
                    det = J[0][0] * J[1][1] * J[2][2] +
                          J[1][0] * J[2][1] * J[0][2] +
                          J[2][0] * J[0][1] * J[1][2] -
                          J[2][0] * J[1][1] * J[0][2] -
                          J[1][0] * J[0][1] * J[2][2] -
                          J[0][0] * J[2][1] * J[1][2];

                    // B and its transpose matrix
                    std::vector<std::vector<Float>> B(6,
                                                      std::vector<Float>(24));
                    std::vector<std::vector<Float>> BT(24,
                                                       std::vector<Float>(6));

                    for (int l = 0; l < 8; l++) {
                        Float dNdx = invJ[0][0] * dNdr1[l] +
                                     invJ[0][1] * dNdr2[l] + invJ[0][2];
                        Float dNdy = invJ[1][0] * dNdr1[l] +
                                     invJ[1][1] * dNdr2[l] + invJ[1][2];
                        Float dNdz = invJ[2][0] * dNdr1[l] +
                                     invJ[2][1] * dNdr2[l] + invJ[2][2];

                        B[0][3 * l] = dNdx;
                        B[1][3 * l] = 0.0;
                        B[2][3 * l] = 0.0;
                        B[3][3 * l] = dNdy;
                        B[4][3 * l] = 0.0;
                        B[5][3 * l] = dNdz;
                        B[0][3 * l + 1] = 0.0;
                        B[1][3 * l + 1] = dNdy;
                        B[2][3 * l + 1] = 0.0;
                        B[3][3 * l + 1] = dNdx;
                        B[4][3 * l + 1] = dNdz;
                        B[5][3 * l + 1] = 0.0;
                        B[0][3 * l + 2] = 0.0;
                        B[1][3 * l + 2] = 0.0;
                        B[2][3 * l + 2] = dNdz;
                        B[3][3 * l + 2] = 0.0;
                        B[4][3 * l + 2] = dNdy;
                        B[5][3 * l + 2] = dNdx;

                        BT[3 * l][0] = dNdx;
                        BT[3 * l][1] = 0.0;
                        BT[3 * l][2] = 0.0;
                        BT[3 * l][3] = dNdy;
                        BT[3 * l][4] = 0.0;
                        BT[3 * l][5] = dNdz;
                        BT[3 * l + 1][0] = 0.0;
                        BT[3 * l + 1][1] = dNdy;
                        BT[3 * l + 1][2] = 0.0;
                        BT[3 * l + 1][3] = dNdx;
                        BT[3 * l + 1][4] = dNdz;
                        BT[3 * l + 1][5] = 0.0;
                        BT[3 * l + 2][0] = 0.0;
                        BT[3 * l + 2][1] = 0.0;
                        BT[3 * l + 2][2] = dNdz;
                        BT[3 * l + 2][3] = 0.0;
                        BT[3 * l + 2][4] = dNdy;
                        BT[3 * l + 2][5] = dNdx;
                    }
                    // stress matrix
                    for (int l = 0; l < 6; l++) {
                        for (int m = 0; m < 24; m++) {
                            tmp[] = B[] * F[];
                        }
                    }
                    for (int l = 0; l < 6; l++) {
                        for (int m = 0; m < 6; m++) {
                            result[] = D[] * tmp[];
                        }
                    }
                }
            }
        }
    }
    return result;
}

std::vector<std::shared_ptr<Shape>> CreateSimpleModelShape(
    const Transform *o2w, const Transform *w2o, bool reverseOrientation,
    const ParamSet &params,
    std::map<std::string, std::shared_ptr<Texture<Float>>> *floatTextures) {
    int nvi, npi, nuvi, nsi, nni;

    Float width = params.FindOneFloat("width", 1.0);
    Float height = params.FindOneFloat("height", 1.0);
    Float depth = params.FindOneFloat("depth", 1.0);
    Float h = params.FindOneFloat("h", 0.1);
    int u = params.FindOneFloat("u", (int)width / h);
    int v = params.FindOneFloat("v", (int)height / h);
    int w = params.FindOneFloat("w", (int)depth / h);
    const int nTriangles = 12;
    const int vertexIndices[3 * nTriangles] = {
        0, 1, 2, 0, 2, 3, 1, 5, 6, 1, 6, 2, 5, 4, 7, 5, 7, 6,
        4, 0, 3, 4, 3, 7, 0, 5, 1, 0, 4, 5, 3, 2, 6, 3, 6, 7};
    const int nVertices = 8;

    const int *vi = params.FindInt("indices", &nvi);
    const Point3f *P = params.FindPoint3f("P", &npi);
    const Point2f *uvs = params.FindPoint2f("uv", &nuvi);
    if (!uvs) uvs = params.FindPoint2f("st", &nuvi);
    std::vector<Point2f> tempUVs;
    if (!uvs) {
        const Float *fuv = params.FindFloat("uv", &nuvi);
        if (!fuv) fuv = params.FindFloat("st", &nuvi);
        if (fuv) {
            nuvi /= 2;
            tempUVs.reserve(nuvi);
            for (int i = 0; i < nuvi; ++i)
                tempUVs.push_back(Point2f(fuv[2 * i], fuv[2 * i + 1]));
            uvs = &tempUVs[0];
        }
    }
    if (uvs) {
        if (nuvi < npi) {
            Error(
                "Not enough of \"uv\"s for triangle mesh.  Expected %d, "
                "found %d.  Discarding.",
                npi, nuvi);
            uvs = nullptr;
        } else if (nuvi > npi)
            Warning(
                "More \"uv\"s provided than will be used for triangle "
                "mesh.  (%d expcted, %d found)",
                npi, nuvi);
    }
    if (!vi) {
        Error(
            "Vertex indices \"indices\" not provided with triangle mesh shape");
        return std::vector<std::shared_ptr<Shape>>();
    }
    if (!P) {
        Error("Vertex positions \"P\" not provided with triangle mesh shape");
        return std::vector<std::shared_ptr<Shape>>();
    }
    const Vector3f *S = params.FindVector3f("S", &nsi);
    if (S && nsi != npi) {
        Error("Number of \"S\"s for triangle mesh must match \"P\"s");
        S = nullptr;
    }
    const Normal3f *N = params.FindNormal3f("N", &nni);
    if (N && nni != npi) {
        Error("Number of \"N\"s for triangle mesh must match \"P\"s");
        N = nullptr;
    }
    for (int i = 0; i < nvi; ++i)
        if (vi[i] >= npi) {
            Error(
                "trianglemesh has out of-bounds vertex index %d (%d \"P\" "
                "values were given",
                vi[i], npi);
            return std::vector<std::shared_ptr<Shape>>();
        }

    int nfi;
    const int *faceIndices = params.FindInt("faceIndices", &nfi);
    if (faceIndices && nfi != nvi / 3) {
        Error("Number of face indices, %d, doesn't match number of faces, %d",
              nfi, nvi / 3);
        faceIndices = nullptr;
    }

    std::shared_ptr<Texture<Float>> alphaTex;
    std::string alphaTexName = params.FindTexture("alpha");
    if (alphaTexName != "") {
        if (floatTextures->find(alphaTexName) != floatTextures->end())
            alphaTex = (*floatTextures)[alphaTexName];
        else
            Error("Couldn't find float texture \"%s\" for \"alpha\" parameter",
                  alphaTexName.c_str());
    } else if (params.FindOneFloat("alpha", 1.f) == 0.f)
        alphaTex.reset(new ConstantTexture<Float>(0.f));

    std::shared_ptr<Texture<Float>> shadowAlphaTex;
    std::string shadowAlphaTexName = params.FindTexture("shadowalpha");
    if (shadowAlphaTexName != "") {
        if (floatTextures->find(shadowAlphaTexName) != floatTextures->end())
            shadowAlphaTex = (*floatTextures)[shadowAlphaTexName];
        else
            Error(
                "Couldn't find float texture \"%s\" for \"shadowalpha\" "
                "parameter",
                shadowAlphaTexName.c_str());
    } else if (params.FindOneFloat("shadowalpha", 1.f) == 0.f)
        shadowAlphaTex.reset(new ConstantTexture<Float>(0.f));

    return CreateTriangleMesh(o2w, w2o, reverseOrientation, nTriangles,
                              vertexIndices, nVertices, P, S, N, uvs, alphaTex,
                              shadowAlphaTex, faceIndices);
}

}  // namespace pbrt
