using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BrukerFileReader
{
    public enum BrukerBafVariableNames
    {
        // Below lines are directly from the bruker baf file sqlite table "SupportedVariables"

        //Digitizer_Summation,2
        //Collision_Energy_Act,5
        //MSMS_PreCursorChargeState,6
        //MSMS_IsolationMass_Act,7
        //Quadrupole_IsolationResolution_Act,8
        //TOF_DeviceTempCurrentValue1,9
        //TOF_DeviceTempCurrentValue2,10
        //TOF_SwitchTempCompensation,11
        //TOF_UseOnlyOneTemperature,12
        //TOF_DeviceTempGradient1,13
        //TOF_DeviceTempGradient2,14
        //TOF_DeviceReferenceTemp1,15
        //TOF_DeviceReferenceTemp2,16
        //Prototype_AutoMSMS_Exclusion_Mass_Count,20
        //Prototype_AutoMSMS_Exclusion_Mass_Bin_Width,21
        //Prototype_AutoMSMS_Exclusion_Mass_List_Blocks,22

        // Names are exactly equivalent to what the Bruker file gives you.
        // Do not change them. 
        Digitizer_Summation = 2,
        Collision_Energy_Act = 5,
        MSMS_PreCursorChargeState = 6,
        MSMS_IsolationMass_Act = 7,
        Quadrupole_IsolationResolution_Act = 8,
        TOF_DeviceTempCurrentValue1 = 9,
        TOF_DeviceTempCurrentValue2 = 10,
        TOF_SwitchTempCompensation = 11,
        TOF_UseOnlyOneTemperature = 12,
        TOF_DeviceTempGradient1 = 13,
        TOF_DeviceTempGradient2 = 14,
        TOF_DeviceReferenceTemp1 = 15,
        TOF_DeviceReferenceTemp2 = 16,
        Prototype_AutoMSMS_Exclusion_Mass_Count = 20,
        Prototype_AutoMSMS_Exclusion_Mass_Bin_Width = 21,
        Prototype_AutoMSMS_Exclusion_Mass_List_Blocks = 22
    }
}
