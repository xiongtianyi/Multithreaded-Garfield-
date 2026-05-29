#include "GPUFunctions.h"
#define __GPUCOMPILE__
#include "ComponentGPU.h"
#undef __GPUCOMPILE__
#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/ComponentElmer.hh"


#define __GPUCOMPILE__

namespace Garfield {

    #include "ComponentFieldMap.cc"
    #include "ComponentAnsys123.cc"

    double Component::CreateGPUTransferObject(ComponentGPU *&comp_gpu)
    {
        // copy periodicity and min/max
        for (int i = 0; i < 3; i++)
        {
            comp_gpu->m_periodic[i] = m_periodic[i];
            comp_gpu->m_mirrorPeriodic[i] = m_mirrorPeriodic[i];
            comp_gpu->m_axiallyPeriodic[i] = m_axiallyPeriodic[i];
            comp_gpu->m_rotationSymmetric[i] = m_rotationSymmetric[i];
        }
        comp_gpu->m_ready = m_ready;

        comp_gpu->m_ComponentType = ComponentGPU::ComponentType::Component;

        return 0;
    }

    double ComponentFieldMap::CreateGPUTransferObject(ComponentGPU *&comp_gpu)
    {
        double alloc{0};

        // Check if bounding boxes of elements have been computed
        if (!m_cacheElemBoundingBoxes) {
            std::cout << m_className << "::CreateGPUTransferObject:\n"
                    << "    Caching the bounding boxes of all elements...";
            CalculateElementBoundingBoxes();
            std::cout << " done.\n";
            m_cacheElemBoundingBoxes = true;
        }

        // initialise the tetraheral tree
        if (m_useTetrahedralTree) {
            if (!m_octree) {
                if (!InitializeTetrahedralTree()) {
                    std::cerr << m_className << "::FindElement13:\n";
                    std::cerr << "    Tetrahedral tree initialization failed.\n";
                    return alloc;
                }
            }
        }

        // copy periodicity and min/max
        for (int i = 0; i < 3; i++)
        {
            comp_gpu->m_mapmin[i] = m_mapmin[i];
            comp_gpu->m_mapmax[i] = m_mapmax[i];
            comp_gpu->m_mapamin[i] = m_mapamin[i];
            comp_gpu->m_mapamax[i] = m_mapamax[i];
        }

        // materials
        comp_gpu->numMaterials = m_materials.size();
        checkCudaErrors(cudaMallocManaged(&(comp_gpu->m_materials), sizeof(ComponentGPU::Material) * comp_gpu->numMaterials));
        alloc += sizeof(ComponentGPU::Material) * comp_gpu->numMaterials;

        for (int i = 0; i < comp_gpu->numMaterials; i++)
        {
            comp_gpu->m_materials[i].eps = m_materials[i].eps;
            comp_gpu->m_materials[i].ohm = m_materials[i].ohm;
            comp_gpu->m_materials[i].driftmedium = m_materials[i].driftmedium;

            if (m_materials[i].medium)
            {
                alloc += m_materials[i].medium->CreateGPUTransferObject(comp_gpu->m_materials[i].medium);
            }
        }

        // tetrahedral tree related things
        comp_gpu->m_checkMultipleElement = m_checkMultipleElement;
        comp_gpu->m_useTetrahedralTree = m_useTetrahedralTree;

        // elements
        comp_gpu->numElements = m_elements.size();
        checkCudaErrors(cudaMallocManaged(&(comp_gpu->m_elements), sizeof(ComponentGPU::Element) * comp_gpu->numElements));
        alloc += sizeof(ComponentGPU::Element) * comp_gpu->numElements;

        checkCudaErrors(cudaMallocManaged(&(comp_gpu->m_degenerate), sizeof(bool) * comp_gpu->numElements));
        alloc += sizeof(bool) * comp_gpu->numElements;

        checkCudaErrors(cudaMallocManaged(&(comp_gpu->m_bbMin), sizeof(GPUFLOAT*) * comp_gpu->numElements));
        checkCudaErrors(cudaMallocManaged(&(comp_gpu->m_bbMax), sizeof(GPUFLOAT*) * comp_gpu->numElements));
        alloc += sizeof(GPUFLOAT*) * comp_gpu->numElements * 2;

        for (int i = 0; i < comp_gpu->numElements; i++)
        {
            checkCudaErrors(cudaMallocManaged(&comp_gpu->m_bbMin[i], sizeof(GPUFLOAT) * 3));
            checkCudaErrors(cudaMallocManaged(&comp_gpu->m_bbMax[i], sizeof(GPUFLOAT) * 3));
            alloc += sizeof(GPUFLOAT) * 3 * 2;
            for (int j = 0; j < 3; j++)
            {
            comp_gpu->m_bbMin[i][j] = m_bbMin[i][j];
            comp_gpu->m_bbMax[i][j] = m_bbMax[i][j];
            }
        }

        for (int i = 0; i < comp_gpu->numElements; i++)
        {
            for (int j = 0; j < 10; j++)
            {
                comp_gpu->m_elements[i].emap[j] = m_elements[i].emap[j];
            }
            comp_gpu->m_elements[i].matmap = m_elements[i].matmap;
            comp_gpu->m_degenerate[i] = m_degenerate[i];
        }

        comp_gpu->numNodes = m_nodes.size();
        checkCudaErrors(cudaMallocManaged(&(comp_gpu->m_nodes), sizeof(ComponentGPU::Node) * comp_gpu->numNodes));
        alloc += sizeof(ComponentGPU::Node) * comp_gpu->numNodes;

        for (int i = 0; i < comp_gpu->numNodes; i++)
        {
            comp_gpu->m_nodes[i].x = m_nodes[i].x;
            comp_gpu->m_nodes[i].y = m_nodes[i].y;
            comp_gpu->m_nodes[i].z = m_nodes[i].z;
        }

        comp_gpu->m_numpot = m_pot.size();
        checkCudaErrors(cudaMallocManaged(&(comp_gpu->m_pot), sizeof(double) * comp_gpu->m_numpot));
        alloc += sizeof(double) * comp_gpu->m_numpot;

        for (int i = 0; i < comp_gpu->m_numpot; i++)
        {
            comp_gpu->m_pot[i] = m_pot[i];
        }

        checkCudaErrors(cudaMallocManaged(&(comp_gpu->m_w12), sizeof(GPUFLOAT**) * comp_gpu->numElements));
        alloc += sizeof(GPUFLOAT**) * comp_gpu->numElements;

        for (int i = 0; i < comp_gpu->numElements; i++)
        {
            checkCudaErrors(cudaMallocManaged(&comp_gpu->m_w12[i], sizeof(GPUFLOAT*) * 4));
            alloc += sizeof(GPUFLOAT*) * 4;

            for (int j = 0; j < 4; j++)
            {
                checkCudaErrors(cudaMallocManaged(&comp_gpu->m_w12[i][j], sizeof(GPUFLOAT) * 3));
                alloc += sizeof(GPUFLOAT) * 3;

                for (int k = 0; k < 3; k++)
                {
                    comp_gpu->m_w12[i][j][k] = m_w12[i][j][k];
                }
            }
        }

        // Weighting potentials
        // - Label information is not used here and instead we are relying on the map being ordered
        comp_gpu->m_num_wpots = m_wpot.size();
        checkCudaErrors(cudaMallocManaged(&(comp_gpu->m_num_entries_wpot), sizeof(int) * comp_gpu->m_num_wpots));
        alloc += sizeof(int) * comp_gpu->m_num_wpots;
        checkCudaErrors(cudaMallocManaged(&(comp_gpu->m_wpot), sizeof(double*) * comp_gpu->m_num_wpots));
        alloc += sizeof(double*) * comp_gpu->m_num_wpots;
        size_t index = 0;
        for (const auto& data: m_wpot) {
            comp_gpu->m_num_entries_wpot[index] = data.second.size();
            checkCudaErrors(cudaMallocManaged(&(comp_gpu->m_wpot[index]), sizeof(double) * data.second.size()));
            alloc += sizeof(double) * data.second.size();
            for (size_t j = 0; j < data.second.size(); j++) {
                comp_gpu->m_wpot[index][j] = data.second.at(j);
            }
            ++index;
        }

        if (comp_gpu->m_useTetrahedralTree)
        {
            alloc += m_octree->CreateGPUTransferObject(comp_gpu->m_octree);
        }

        comp_gpu->numElements = m_elements.size();
        // TODO TN GPU: Fix type to CurvedTetrahedron, work out how to get this
        // from the CPU FieldMap in future
        comp_gpu->m_elementType = ComponentGPU::ElementType::CurvedTetrahedron;
        comp_gpu->m_is3d = m_is3d;

        comp_gpu->m_ComponentType = ComponentGPU::ComponentType::ComponentFieldMap;

        return alloc;
    }

    double ComponentAnsys123::CreateGPUTransferObject(ComponentGPU *&comp_gpu)
    {
        // create main sensor GPU class
        checkCudaErrors(cudaMallocManaged(&comp_gpu, sizeof(ComponentGPU)));
        double alloc{sizeof(ComponentGPU)};

        alloc += ComponentFieldMap::CreateGPUTransferObject(comp_gpu);
        alloc += Component::CreateGPUTransferObject(comp_gpu);

        comp_gpu->m_ComponentType = ComponentGPU::ComponentType::ComponentAnsys123;
        return alloc;
    }

    double ComponentElmer::CreateGPUTransferObject(ComponentGPU *&comp_gpu)
    {
        // create main sensor GPU class
        checkCudaErrors(cudaMallocManaged(&comp_gpu, sizeof(ComponentGPU)));
        double alloc{sizeof(ComponentGPU)};

        alloc += ComponentFieldMap::CreateGPUTransferObject(comp_gpu);
        alloc += Component::CreateGPUTransferObject(comp_gpu);

        comp_gpu->m_ComponentType = ComponentGPU::ComponentType::ComponentElmer;
        return alloc;
    }

    double ComponentComsol::CreateGPUTransferObject(ComponentGPU *&comp_gpu)
    {
        // create main sensor GPU class
        checkCudaErrors(cudaMallocManaged(&comp_gpu, sizeof(ComponentGPU)));
        double alloc{sizeof(ComponentGPU)};

        alloc += ComponentFieldMap::CreateGPUTransferObject(comp_gpu);
        alloc += Component::CreateGPUTransferObject(comp_gpu);

        comp_gpu->m_ComponentType = ComponentGPU::ComponentType::ComponentComsol;
        return alloc;
    }
}
