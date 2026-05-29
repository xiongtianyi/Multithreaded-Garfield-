#define __GPUCOMPILE__
#include "TetrahedralTree.cc"
#undef __GPUCOMPILE__
#include "Garfield/TetrahedralTree.hh"

namespace Garfield {

    double TetrahedralTree::CreateGPUTransferObject(TetrahedralTreeGPU *&tree_gpu)
    {
        // create main TetrahedralTree GPU class
        checkCudaErrors(cudaMallocManaged(&tree_gpu, sizeof(TetrahedralTreeGPU)));
        double alloc{sizeof(TetrahedralTreeGPU)};

        // copy the elements
        tree_gpu->numelements = elements.size();
        if (elements.size() != 0)
        {
            alloc += CreateGPUArrayFromVector<int>(elements, tree_gpu->numelements, tree_gpu->elements);
        }
        else
        {
            tree_gpu->elements = nullptr;
        }

        // copy children
        for (int i = 0; i < 8; i++)
        {
            if (children[i])
            {
                alloc += children[i]->CreateGPUTransferObject(tree_gpu->children[i]);
            }
            else
            {
                tree_gpu->children[i] = nullptr;
            }
        }

        // copy other vars
        tree_gpu->m_origin = Vec3GPU{m_origin.x, m_origin.y, m_origin.z};
        tree_gpu->m_halfDimension = Vec3GPU{m_halfDimension.x, m_halfDimension.y, m_halfDimension.z};
        tree_gpu->m_min = Vec3GPU{m_min.x, m_min.y, m_min.z};
        tree_gpu->m_max = Vec3GPU{m_max.x, m_max.y, m_max.z};

        return alloc;
    }
}