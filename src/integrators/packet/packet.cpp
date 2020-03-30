/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/ray_sse.h>
#include <array>

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength("Packet tracer", "Average path length", EAverage);


class MIPacketTracer : public MCPacketIntegrator {
public:
    MIPacketTracer(const Properties &props)
        : MCPacketIntegrator(props) { }

    /// Unserialize from a binary data stream
    MIPacketTracer(Stream *stream, InstanceManager *manager)
        : MCPacketIntegrator(stream, manager) { }

    std::array<Spectrum, 4> Li_packet(RayPacket4 &packet, RadianceQueryRecord *rRec) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec[0].scene;
        Intersection &its = rRec[0].its;
        RayDifferential ray;
        std::array<Spectrum, 4> Lis;

        RayInterval4 rayIntervals;
        Intersection4 itss;

        bool scattered = false;

        const RayPacket4 packet_ = packet;

        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        //rRec.rayIntersect(ray);

        scene->rayIntersectPacketIncoherent(packet_, rayIntervals, itss, nullptr);

        //ray.mint = Epsilon;

        

        /* Store statistics */
        avgPathLength.incrementBase();
        avgPathLength += rRec[1].depth;

        return std::array<Spectrum, 4>();
    }

    inline Float miWeight(Float pdfA, Float pdfB) const {
        pdfA *= pdfA;
        pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        MCPacketIntegrator::serialize(stream, manager);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "MIPacketTracer[" << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  strictNormals = " << m_strictNormals << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(MIPacketTracer, false, MCPacketIntegrator)
MTS_EXPORT_PLUGIN(MIPacketTracer, "MI packet tracer");
MTS_NAMESPACE_END
