#include "GPUFunctions.h"
#define __GPUCOMPILE__
#include "MediumGPU.h"
#undef __GPUCOMPILE__
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/GarfieldConstants.hh"
#include "RandomGPU.h"

#define __GPUCOMPILE__

namespace Garfield {

    #include "MediumMagboltz.cc"

    __device__ double MediumGPU::GetElectronCollisionRate(const GPUFLOAT e,
        const int band)
    {
        switch(m_MediumType)
        {
        case MediumType::Medium:
            {
                return 0;
            }
        case MediumType::MediumMagboltz:
            {
                return GetElectronCollisionRate__MediumMagboltz(e, band);
            }

        }

        return 0;
    }

    double Medium::CreateGPUTransferObject(MediumGPU *&med_gpu)
    {
        double alloc{0};

        med_gpu->m_driftable = m_driftable;
        med_gpu->m_ionisable = m_ionisable;
        med_gpu->m_id = m_id;
        med_gpu->m_microscopic = m_microscopic;

        med_gpu->m_MediumType = MediumGPU::MediumType::Medium;
        return alloc;
    }

    double MediumGas::CreateGPUTransferObject(MediumGPU *&med_gpu)
    {
        double alloc{0};
        med_gpu->m_MediumType = MediumGPU::MediumType::MediumGas;
        return alloc;
    }

    double MediumMagboltz::CreateGPUTransferObject(MediumGPU *&med_gpu)
    {
        checkCudaErrors(cudaMallocManaged(&med_gpu, sizeof(MediumGPU)));
        double alloc{sizeof(MediumGPU)};

        alloc += Medium::CreateGPUTransferObject(med_gpu);
        alloc += MediumGas::CreateGPUTransferObject(med_gpu);

        med_gpu->m_eHighLog = m_eHighLog;
        med_gpu->m_lnStep = m_lnStep;
        med_gpu->m_eHigh = m_eHigh;
        med_gpu->m_eStep = m_eStep;
        med_gpu->m_eStepInv = m_eStepInv;
        med_gpu->m_nTerms = m_nTerms;
        med_gpu->m_useOpalBeaty = m_useOpalBeaty;
        med_gpu->m_useAnisotropic = m_useAnisotropic;
        med_gpu->m_cfNull = m_cfNull;
        med_gpu->m_eMax = m_eMax;

        alloc += CreateGPUArrayOfArraysFromVector<double>(m_cf, med_gpu->m_cf, med_gpu->m_numcf, med_gpu->m_numcfIdx);
        alloc += CreateGPUArrayOfArraysFromVector<double>(m_cfLog, med_gpu->m_cfLog, med_gpu->m_numcfLog, med_gpu->m_numcfLogIdx);
        alloc += CreateGPUArrayOfArraysFromVector<double>(m_scatPar, med_gpu->m_scatPar, med_gpu->m_numscatPar, med_gpu->m_numscatParIdx);
        alloc += CreateGPUArrayOfArraysFromVector<double>(m_scatParLog, med_gpu->m_scatParLog, med_gpu->m_numscatParLog, med_gpu->m_numscatParLogIdx);
        alloc += CreateGPUArrayOfArraysFromVector<double>(m_scatCut, med_gpu->m_scatCut, med_gpu->m_numscatCut, med_gpu->m_numscatCutIdx);
        alloc += CreateGPUArrayOfArraysFromVector<double>(m_scatCutLog, med_gpu->m_scatCutLog, med_gpu->m_numscatCutLog, med_gpu->m_numscatCutLogIdx);

        for (int i = 0; i < Magboltz::nMaxLevels; i++)
        {
            med_gpu->m_csType[i] = m_csType[i];
            med_gpu->m_wOpalBeaty[i] = m_wOpalBeaty[i];
            med_gpu->m_yFluorescence[i] = m_yFluorescence[i];
            med_gpu->m_nAuger1[i] = m_nAuger1[i];
            med_gpu->m_nAuger2[i] = m_nAuger2[i];
            med_gpu->m_eAuger1[i] = m_eAuger1[i];
            med_gpu->m_eAuger2[i] = m_eAuger2[i];
            med_gpu->m_nFluorescence[i] = m_nFluorescence[i];
            med_gpu->m_eFluorescence[i] = m_eFluorescence[i];
            med_gpu->m_scatModel[i] = m_scatModel[i];
            med_gpu->m_energyLoss[i] = m_energyLoss[i];
        }

        for (int i = 0; i < m_nMaxGases; i++)
        {
            med_gpu->m_rgas[i] = m_rgas[i];
            med_gpu->m_s2[i] = m_s2[i];
        }

        alloc += CreateGPUArrayFromVector<double>(m_cfTot, med_gpu->m_numcfTot, med_gpu->m_cfTot);
        alloc += CreateGPUArrayFromVector<double>(m_cfTotLog, med_gpu->m_numcfTotLog, med_gpu->m_cfTotLog);

        med_gpu->m_MediumType = MediumGPU::MediumType::MediumMagboltz;

        return alloc;
    }



}