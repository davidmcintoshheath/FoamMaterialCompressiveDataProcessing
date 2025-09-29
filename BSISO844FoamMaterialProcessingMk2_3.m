clc;
clear all;
%find all the data files
tic
fileNames = findDirectoryTextFileNames();
vectorOfMaterialData = [];
vectorOfDataForAveraging = zeros(1,3*length(fileNames));
moveToNextDataStream = 0;
shortestVector = -1;
writeOutPutFile = "Mk2_3ABAQUSSawbonesMaterialDataTest";
writeOutPutFile = createDataFileName(fileNames(1),writeOutPutFile);
for fileNameIterator = 1:length(fileNames)
    %extrct the file name data
    nameData = fileNameDataExtraction(fileNames(fileNameIterator));
    foamDimentionData = [];
    foamDimentionCounter = 1;
    for(dimentionIterator = 1:2:(length(nameData)-1))
        foamDimentionData(foamDimentionCounter) = (nameData(dimentionIterator) + nameData(dimentionIterator+1))/2;
        foamDimentionCounter = foamDimentionCounter + 1;
    end
    %extract the file data
    foamCompressionData = [];
    foamCompressionData = dataTextFileToRawData(fileNames(fileNameIterator));
    foamDataSIConvertDataStreamArray = [2];
    foamCompressionData = convertToSIUnits(foamCompressionData,foamDataSIConvertDataStreamArray);
    foamDimentionDataSelectionArray = [1,2,3];
    foamDimentionData = convertToSIUnits(foamDimentionData,foamDimentionDataSelectionArray);
    foamCompressionData = convertForceDisplacementToStressStrainData(foamCompressionData,foamDimentionData);
    %process the data to work out the upper stress level and calculate the modulus
    processedModulusDataStruct = generateStructOfKeyFoamDataPoints(foamCompressionData); %Stress strain in engineering stress strain
    vectorOfMaterialData(fileNameIterator,1) = processedModulusDataStruct.modulus;
    vectorOfMaterialData(fileNameIterator,2) = processedModulusDataStruct.ultimateCompressiveStress;
    vectorOfMaterialData(fileNameIterator,3) = processedModulusDataStruct.compressiveYieldStress;
    
    %Correct the data for the linear region
    [foamCompressionData,processedModulusDataStruct] = correctLinearRegion(foamCompressionData, processedModulusDataStruct);
    
    %This converts the foam engineering stress data to true stress strain
    %data. Has an error of ~0.024% with relation to converting the final
    %average stress strain graph to true and then doing the calculations
    %Move the true foam stress calculation to after the averaging step
    %foamTrueStressStrainData = convertEngineeringStressDataToTrueStressData(foamCompressionData);
    foamTrueStressStrainData = foamCompressionData;
    
    %Add to averaging Vector
    for(addToAveragingVectorIterator = 1:length(foamTrueStressStrainData))
        for(dataStreamSelector = 1:3)
            vectorOfDataForAveraging(addToAveragingVectorIterator,(((dataStreamSelector-1)*length(fileNames) +1)+moveToNextDataStream)) = foamTrueStressStrainData(addToAveragingVectorIterator, dataStreamSelector);
        end
    end
    moveToNextDataStream = moveToNextDataStream + 1;
    
    %hold on
    %plot(foamCompressionData(:,2),foamCompressionData(:,3),'b');
    %plot(trueFoamCompressionData(:,2),trueFoamCompressionData(:,3));
    %plot(plasticStrainData(:,2),plasticStrainData(:,3));
    %plot(trueFoamCompressionData(:,2),trueFoamCompressionData(:,3));

    %Consolidation phase is the least critical for the success of the
    %simulation so could deal with the most error. Cut the tails off the
    %data

    %Convert to plastic (Changes UCS by, where as the otherwise 0.04%
    %change by averaging)

    %Work out the new "TRUE" ucs and yeild
    %Convert to plastic (Changes UCS by, where as the otherwise 0.04%
    %change by averaging)
    
    %Generate a series of data points including the yeild UCS and End
    %points Distribute in (10/20 before UCS, 6/20 in 0.75 msx strain, 4/20 last 0.25
    %max strain) 
    % This actually caused a major issue in ABAQUS as ABAQUS wants to regularize the 
    % data to create a table with equal strain steps. This is apparently done for 
    % performance. Becue the strain steps created by this algorithm
    % increase by orders of magnitude from start to end of the material curve, ABAQUS 
    % throws a "Regularization Error". 

    %Strain was converted after averaging and has had a new property struct
    %generated. This increased the modulus by about 10% but represents the
    %true material data better than the engineering stress strain.

    %What is left to do:  
    
    %Fix the data averageing and fix the data stream length equalization
    %functions
end

%average the strain and stress data over the material curves
strainAveragingStepSize = 0.0008; %was 0.004
averagedMaterialCurveData = stressStrainCurveAverage(vectorOfDataForAveraging,strainAveragingStepSize);

%average the material properties of the data curves
averageVectorOfMaterialData = averageMaterialProperties(vectorOfMaterialData);

%correct the material data struct for properties
averagedMaterialCurveData = convertEngineeringStressDataToTrueStressData(averagedMaterialCurveData);

%Decided to run a curve on the true stress strain because the true stress strain conversion does change the modulus
% and it seems reasnable that the material model using the true stress and strain should reflect that
unifiedAverageMaterialCurveStruct  = generateStructOfKeyFoamDataPoints(averagedMaterialCurveData);


%Correct the linear Region
[averagedMaterialCurveData,unifiedAverageMaterialCurveStruct] = correctLinearRegion(averagedMaterialCurveData, unifiedAverageMaterialCurveStruct);

%plot the new curves
hold on
plot(averagedMaterialCurveData(:,2),averagedMaterialCurveData(:,3),'b');

%output Plastic strain data
plasticStrainData = convertStrainToAbaqusPlasticStrain(averagedMaterialCurveData, unifiedAverageMaterialCurveStruct);
plot(plasticStrainData(:,2),plasticStrainData(:,3),'r');


%Select out the data for input to ABAQUS
numberOfDataPoints = 100;
plasticDataToBeWritten = extractEquallyDistributedABAQUSPlasticStrainDataPoints(plasticStrainData, unifiedAverageMaterialCurveStruct, numberOfDataPoints);
plot(plasticDataToBeWritten(:,1),plasticDataToBeWritten(:,2));

%Intergrate to find the fracture strain energy
[unifiedAverageMaterialCurveStruct] = getABAQUSPlasticFractureEnergy(plasticDataToBeWritten, unifiedAverageMaterialCurveStruct);

%Write data out to text file in both pythonVersion and transfurable
%version
writeABAQUSPythonDataToFile(writeOutPutFile, plasticDataToBeWritten, unifiedAverageMaterialCurveStruct);
writeABAQUSDataToFile(writeOutPutFile, plasticDataToBeWritten, unifiedAverageMaterialCurveStruct);



%place holder assignment
toc
EndOfProgramRunHold = 1;

%%Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data extraction functions
%Data about files in location
function [fileNames] = findDirectoryTextFileNames()
    files = dir("*.txt*");
    fileNames = string({files.name});
end
%Data from file name
function [extractedFileNameData] = fileNameDataExtraction(fileName)
    extractedFileNameData = [];
    endOfFileNameDataStream = 1;
    assignmentCounter = 1;
    stringData = char(fileName);
    searchCharacterVector = [];
    searchCharacterVector = strfind(fileName,' ');
    for(searchIterator = 1:(length(searchCharacterVector)-1))
        extractedCharacterVector = 'ab';
        for(characterAssignmentIterator = (searchCharacterVector(searchIterator)+1):(searchCharacterVector(searchIterator+1)-1))
            extractedCharacterVector(characterAssignmentIterator - searchCharacterVector(searchIterator)) = stringData(characterAssignmentIterator);
        end
        extractedFileNameData(assignmentCounter) = str2double(extractedCharacterVector);
        assignmentCounter = assignmentCounter + 1;
    end
    sizeOfExtractedData = size(extractedFileNameData);
    assert(sizeOfExtractedData(1)~=0,("[DATA ERROR]: THERE WAS NO DATA IN THE FILE NAME"));
end

%Data from text file
function [fileDataOut] = dataTextFileToRawData(fileName)
    fileID = fopen(fileName,"r");
    switch(fileID)
        case(-1)
            fprintf(strcat("File: ",fileName," ", "opening failed\n"));
            return;
        case(fileID)
            fprintf(strcat("File: ",fileName," ", "opened successfully\n"));
    end
    temporaryLineData = 0;
    fileDataOut = [];
    rowIteratorVariable = 1;
    fileOutputData = 0;
    checkWritingData = 0;
    while(~feof(fileID))
        endOfLineCheckVariable = 1;
        columnIteratorVariable = 1;
        fileLineStringData = fgetl(fileID);
        while(endOfLineCheckVariable)
            [fileOutputData, fileLineStringData] = strtok(fileLineStringData);
            fileOutputData = str2double(fileOutputData);
            endOfLineCheckVariable = isfinite(fileOutputData);
            switch(endOfLineCheckVariable)
                case(1)
                    checkWritingData = 1;
                    fileDataOut(rowIteratorVariable,columnIteratorVariable) = fileOutputData;
                    columnIteratorVariable = columnIteratorVariable + 1;
                case(0)
                    continue;
            end
        end
        switch(checkWritingData)
            case(1)
                rowIteratorVariable = rowIteratorVariable + 1;
            case(0)
                continue;
        end
    end
    fileCloseComparator = fclose(fileID);
    switch(fileCloseComparator)
        case(-1)
            fprintf(strcat("File: ",fileName," ", "closing failed\n"));
            return;
        case(0)
            fprintf(strcat("File: ",fileName," ", "closed successfully\n"));
    end
end

function [outputString] = createDataFileName(inputCheckString,inputFileNameToBeModified)
    if(strfind(inputCheckString,"PCF12")>0)
        outputString = strcat(inputFileNameToBeModified,"PCF12");
    elseif(strfind(inputCheckString,"PCF15")>0)
        outputString = strcat(inputFileNameToBeModified,"PCF15");
    elseif(strfind(inputCheckString,"PCF20")>0)
        outputString = strcat(inputFileNameToBeModified,"PCF20");
    elseif(strfind(inputCheckString,"PCF25")>0)
        outputString = strcat(inputFileNameToBeModified,"PCF25");
    else
        assert((1==0),"[INPUT FILE NAME ERROR]: THE FILE NAME DOES NOT CONTAIN THE FOAM DENSITY");
    end
    if(strfind(inputCheckString,"1mmpmin")>0)
        outputString = strcat(outputString,"-1mmpmin");
    elseif(strfind(inputCheckString,"10mmpmin")>0)
        outputString = strcat(outputString,"-10mmpmin");
    elseif(strfind(inputCheckString,"100mmpmin")>0)
        outputString = strcat(outputString,"-100mmpmin");
    else
        assert((1==0),"[INPUT FILE NAME ERROR]: THE FILE NAME DOES NOT CONTAIN THE STRAIN RATE");
    end
end

%%Data Processing funcitons%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data selection with direction
function  [parsed_data] = data_selection_with_location(input_data, direction_charecter, data_location)
    assert(~(data_location >= length(input_data)),"[DATA SIZE ERROR]: There is no room for data processing, value too near upper boundary conditions")
    assert(~((data_location == 1) && (direction_charecter == '<')),"[DATA SIZE ERROR]: There is no room for data processing, value too near lower boundary conditions");
    
    parsed_data = [];
    data_size = size(input_data);
    parsed_itterator_value = 1;
    if(direction_charecter == '>')
        for(itterator_value = data_location:length(input_data))
            for(data_itterator_value = 1:data_size(2))
                parsed_data(parsed_itterator_value, data_itterator_value) = input_data(itterator_value, data_itterator_value);
            end
            parsed_itterator_value = parsed_itterator_value + 1;
        end
    else if(direction_charecter == '<')
            for(itterator_value = 1:data_location)
            for(data_itterator_value = 1:data_size(2))
                parsed_data(parsed_itterator_value, data_itterator_value) = input_data(itterator_value, data_itterator_value);
            end
            parsed_itterator_value = parsed_itterator_value + 1;
        end
        end
    end
end

%Data cleanup and invert direction
function [clean_data] = dartec_text_data_cleanup_and_invert(Raw_data,test_speed,time_offset,displacement_offset,load_offset)
    number_of_data_streams = size(Raw_data);
    time_step = (Raw_data(3,1)-Raw_data(1,1))/3;
    cutoff_boundary = 0.011 + (time_step * test_speed); %0.011 value gained from maximum displacement variation betwene two points over 15138 datapoints in a dwell on the DARTEC at 250Hz.
    if(number_of_data_streams(2) == 4)
        clean_data = [];
        comparrison_data_point = 0;
        assignemnt_varariable = 1;
        for(i = 1:length(Raw_data)-2)
            comp_data_outer = -(Raw_data(i+2,2)-Raw_data(1,2));
            if((comparrison_data_point >= (comp_data_outer-cutoff_boundary))&&(comparrison_data_point <= (comp_data_outer+cutoff_boundary)))
                clean_data(assignemnt_varariable,1) = (Raw_data(i+2,1)-time_offset);
                clean_data(assignemnt_varariable,2) = -(Raw_data(i+2,2)-displacement_offset);
                clean_data(assignemnt_varariable,3) = -(Raw_data(i+2,4)-load_offset);
                comparrison_data_point = abs(Raw_data(i+2,2)-Raw_data(1,2));
                assignemnt_varariable = assignemnt_varariable + 1;
            end
        end
        if(isempty(clean_data))
            fprintf(strcat("Data clean up failed"," ","No data found\n"));
        else
            fprintf("Data clean up succsessful\n");
        end
    end
    if(number_of_data_streams(2) == 3)
        clean_data = [];
        comparrison_data_point = 0;
        assignemnt_varariable = 1;
        for(i = 1:length(Raw_data)-2)
            comp_data_outer = -(Raw_data(i+2,2)-Raw_data(1,2));
            if((comparrison_data_point >= (comp_data_outer-cutoff_boundary))&&(comparrison_data_point <= (comp_data_outer+cutoff_boundary)))
                clean_data(assignemnt_varariable,1) = (Raw_data(i+2,1)-time_offset);
                clean_data(assignemnt_varariable,2) = -(Raw_data(i+2,2)-displacement_offset);
                clean_data(assignemnt_varariable,3) = -(Raw_data(i+2,3)-load_offset);
                comparrison_data_point = abs(Raw_data(i+2,2)-Raw_data(1,2));
                assignemnt_varariable = assignemnt_varariable + 1;
            end
        end
        if(isempty(clean_data))
            fprintf(strcat("Data clean up failed"," ","No data found\n"));
        else
            fprintf("Data clean up succsessful\n");
        end
    end
end

%Data translation function
function [translated_data] = data_translating_function(input_data, data_stream_select, translation_value)
    translated_data = [];
    input_data_size = size(input_data);
    for(itterator_variable = 1:length(input_data))
        for(data_stream_itterator_variable = 1:input_data_size(2))
            if(data_stream_itterator_variable == data_stream_select)
                translated_data(itterator_variable,data_stream_itterator_variable) = input_data(itterator_variable,data_stream_itterator_variable) + translation_value;
            else
                translated_data(itterator_variable,data_stream_itterator_variable) = input_data(itterator_variable,data_stream_itterator_variable);
            end
        end
    end
end

%Data trimming function
function [trimmed_data] = data_trimming_function(input_data, data_stream_select, direction_define_string, limit)
    trimmed_data = [];
    trimmed_itterator_variable = 1;
    input_dimentions = size(input_data);
    if(direction_define_string == '>')
        for(itterator_variable = 1:length(input_data))
            if(input_data(itterator_variable, data_stream_select) >= limit)
                for(data_stream_itterator_variable = 1:input_dimentions(2))
                    trimmed_data(trimmed_itterator_variable,data_stream_itterator_variable) = input_data(itterator_variable,data_stream_itterator_variable);
                end
                trimmed_itterator_variable = trimmed_itterator_variable + 1;
            end
        end
    else if(direction_define_string == '<')
        for(itterator_variable = 1:length(input_data))
            if(input_data(itterator_variable, data_stream_select) <= limit)
                for(data_stream_itterator_variable = 1:input_dimentions(2))
                    trimmed_data(trimmed_itterator_variable,data_stream_itterator_variable) = input_data(itterator_variable,data_stream_itterator_variable);
                end
                trimmed_itterator_variable = trimmed_itterator_variable + 1;
            end
        end
    end
    end
end

%Force displacement data to stress strain data conversion function
function [stressStrainData] = convertForceDisplacementToStressStrainData(inputData,foamDimentionalData)
    %Input should be in SI units
    %Input data should have Time, Displacement, Force arrangement only
    %Foam dimentional data should have X Dimention, Y Dimention, Z
    %Dimention arrangement
    %Check Inputs
    inputDataSize = size(inputData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 3), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME, DISPLACEMENT, FORCE");
    %Run calculation
    foamArea = foamDimentionalData(1) * foamDimentionalData(2);
    stressStrainData = zeros(length(inputData),2);
    for(streamIterator = 1:length(inputData))
        stressStrainData(streamIterator,1) = inputData(streamIterator,1);
        stressStrainData(streamIterator,2) = inputData(streamIterator,2)/foamDimentionalData(3);
        stressStrainData(streamIterator,3) = inputData(streamIterator,3)/foamArea;
    end
end

%Convert to SI Units
function [outputDataInSIUnits] = convertToSIUnits(inputData,transformationSelectionArray)
    %Assume input in mm will convert to m
    outputDataInSIUnits = inputData;
    sizeInputData = size(inputData);
    for(conversionIterator = 1:sizeInputData(1))
        for(streamSelectionIterator = 1:length(transformationSelectionArray))
            outputDataInSIUnits(conversionIterator,transformationSelectionArray(streamSelectionIterator)) = inputData(conversionIterator,transformationSelectionArray(streamSelectionIterator))/1000; 
        end
    end
end

%Foam modulus determination points function
function [outputDataStruct] = generateStructOfKeyFoamDataPoints(inputFoamData)
    DEBUG = false;
    %foam data should be input as TIME, STRAIN, STRESS
    %Check Inputs
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 3), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME, DISPLACEMENT, FORCE");
    assert((inputDataSize(1) > 1), "[DATA INPUT ERROR]: TOO FEW DATA POINTS INPUT");
    %Calculation
    
    %Search for the first maximum stress point (UCS)
    ultimateCompressiveStressDataStruct = findUltimateCompressiveStressStruct(inputFoamData);

    %Search for the section
    ultimateCompressionStressIndex = ultimateCompressiveStressDataStruct.ultimateCompressiveStressValueLocation;
    ultimateCompressionStress = ultimateCompressiveStressDataStruct.ultimateCompressiveStressValue;
    upperStressLimit = ultimateCompressionStress*0.75;
    lowerStressLimit = ultimateCompressionStress*0.25;
    lowerLimitAssignmentFlag = 0;
    upperLimitAssignmentFlag = 0;
    lowerLimitStressIndex = 0;
    upperLimitStressIndex = 0;
    for(searchSelectionIterator = 1:length(inputFoamData))
        %SPECIAL CASE SWITCHES INTEND TO HAVE MOST GO OUTSIDE RANGE, DONT
        %NEED OTHERWISE
        switch((inputFoamData(searchSelectionIterator,3) > lowerStressLimit) && (lowerLimitAssignmentFlag == 0))
            case 1
                lowerLimitStressIndex = searchSelectionIterator;
                lowerLimitAssignmentFlag = 1;
        end
        switch((inputFoamData(searchSelectionIterator,3) > upperStressLimit) && (upperLimitAssignmentFlag == 0))
            case 1
                upperLimitStressIndex = searchSelectionIterator;
                upperLimitAssignmentFlag = 1;
        end
        switch((lowerLimitAssignmentFlag == 1) && (upperLimitAssignmentFlag == 1))
            case 1
                break;
        end
    end
    
    %Assign each of the modulus case
    numberOfDevisions = 6;
    numberOfStepsInEachSection = ceil((upperLimitStressIndex-lowerLimitStressIndex) / numberOfDevisions);
    overSelectionSteps = ((numberOfStepsInEachSection * numberOfDevisions) - (upperLimitStressIndex - lowerLimitStressIndex));
    
    %Account for over selection
    lowerLimitStressIndex = lowerLimitStressIndex - ceil(overSelectionSteps / 2);
    upperLimitStressIndex = upperLimitStressIndex + floor(overSelectionSteps / 2);
    
    %Assign the variables
    modulusSelection = zeros(numberOfStepsInEachSection,2,numberOfDevisions);
    assignmentCounter = lowerLimitStressIndex;
    for(vectorIterator = 1:numberOfDevisions)
        for(assignmentIterator = 1:numberOfStepsInEachSection)
            for(streamIterator = 1:2)
                modulusSelection(assignmentIterator,streamIterator,vectorIterator) = inputFoamData(assignmentCounter,streamIterator+1);
            end
            assignmentCounter = assignmentCounter + 1;
        end
    end
    %Create the fits
    modulusCoefficientVector = zeros(numberOfDevisions,2); %Output modulus, +c
    fitTypeDesignation = 'linear';
    for(fittingIterator = 1:numberOfDevisions)
        linearModel = fitlm(modulusSelection(:,1,fittingIterator),modulusSelection(:,2,fittingIterator),fitTypeDesignation);
        modulusCoefficientVector(fittingIterator,1) = table2array(linearModel.Coefficients(2,1));%fitOutput.a;
        modulusCoefficientVector(fittingIterator,2) = table2array(linearModel.Coefficients(1,1));
    end
    
    %Select the maximum modulus
    maximumModulusComparitorVariable = 0;
    maximumModulusComparitorIndex = 0;
    
    for(maximumModulusIterator = 1:(length(modulusCoefficientVector)-1))
        nextAddition = 0;
        for(modulusAddingIterator = 0:1)
            nextAddition = nextAddition + modulusCoefficientVector(maximumModulusIterator+modulusAddingIterator);
        end
        switch(maximumModulusComparitorVariable < nextAddition)
            case 1
                maximumModulusComparitorVariable = nextAddition;
                maximumModulusComparitorIndex = maximumModulusIterator;
            case 0
                continue;
        end
    end
    
    %Use maximum modulus to work out the inflexion point
    modulusCalculationVector = zeros(numberOfStepsInEachSection*2,2);
    dataAssignmentModulusCounterIndex = 1;
    for(modulusCalculationIterator = maximumModulusComparitorIndex:(maximumModulusComparitorIndex+1))
        for(dataAssignmentModulusCalculationIterator = 1:numberOfStepsInEachSection)
            for(assignmentStreamSelectionIterator = 1:2)
                modulusCalculationVector(dataAssignmentModulusCounterIndex, assignmentStreamSelectionIterator) = modulusSelection(dataAssignmentModulusCalculationIterator,assignmentStreamSelectionIterator,modulusCalculationIterator);
            end
            dataAssignmentModulusCounterIndex = dataAssignmentModulusCounterIndex + 1;
        end
    end
    %Work out the failure stress index of the foam
    failureStressIndex = lowerLimitStressIndex + ((maximumModulusComparitorIndex + 1) * numberOfStepsInEachSection) + 1;

    %Assign the linear model parameters to the outputStruct
    outputDataStruct = struct('modulus',0,'modulusIntercept',0,'ultimateCompressiveStress',0,'indexLocationUCS',0,'compressiveYieldStress',0,'compressiveYieldStressIndex',0);
    
    linearModel = fitlm(modulusCalculationVector(:,1),modulusCalculationVector(:,2),fitTypeDesignation);
    
    %Work out the strain axis intercept (Zero stress strain offset)
    %strainOffsetValue = (-table2array(linearModel.Coefficients(1,1))/table2array(linearModel.Coefficients(2,1))); %old way
    strainOffsetValue = inputFoamData(failureStressIndex,2) - (inputFoamData(failureStressIndex,3) / table2array(linearModel.Coefficients(2,1)));
    
    %Find the the closest strain value to the strain offset
    closestStrainValue = 0;
    for(strainOffsetSearchIterator = 2:length(inputFoamData))
        switch( inputFoamData(strainOffsetSearchIterator,2) >= strainOffsetValue )
            case true
                strainStep = inputFoamData(strainOffsetSearchIterator,2) - inputFoamData(strainOffsetSearchIterator-1,2);
                strainDifferenceWithSearchValue = strainOffsetValue - inputFoamData(strainOffsetSearchIterator-1,2);
                if((strainStep / 2) <= strainDifferenceWithSearchValue)
                    closestStrainValue = inputFoamData(strainOffsetSearchIterator,2);
                else
                    closestStrainValue = inputFoamData(strainOffsetSearchIterator - 1,2);
                end
                break;
            case false
                continue;
        end
    end
    strainOffsetIndex = findValueIndex(closestStrainValue,inputFoamData(:,2),1);
    
    %Write data to foam output data struct
    outputDataStruct.modulus = table2array(linearModel.Coefficients(2,1));
    outputDataStruct.modulusIntercept = table2array(linearModel.Coefficients(1,1));
    outputDataStruct.ultimateCompressiveStress = ultimateCompressionStress;
    outputDataStruct.indexLocationUCS = ultimateCompressionStressIndex;
    outputDataStruct.compressiveYieldStress = inputFoamData(failureStressIndex,3);
    outputDataStruct.compressiveYieldStressIndex = failureStressIndex;
    outputDataStruct.strainAtYield = outputDataStruct.compressiveYieldStress / outputDataStruct.modulus;
    outputDataStruct.strainOffsetValue = strainOffsetValue;
    outputDataStruct.strainOffsetClosestIndex= strainOffsetIndex;
    outputDataStruct.plasticToughnessEnergy = 0;
    
    %check output
    logDisp("Checking the accuracy of the model",-1,DEBUG);
    if(outputDataStruct.ultimateCompressiveStress == inputFoamData(outputDataStruct.indexLocationUCS,3))
        logDisp("The ultimate compressive stress index is correct",1,DEBUG);
    else
        logErr("The ultimate compressive stress index is wrong, check",-1,DEBUG);
    end
end
%

%Convert engineering stress strain to true stress strain
function [trueStressStrainData] = convertEngineeringStressDataToTrueStressData(inputFoamData)
    %foam data should be input as TIME, STRAIN, STRESS
    %Check Inputs
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 3), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME, STRAIN, STRESS");
    %Calculation
    trueStressStrainData = zeros(inputDataSize);
    for(stressStrainIterator = 1:inputDataSize(1))
        trueStressStrainData(stressStrainIterator,1) = inputFoamData(stressStrainIterator,1);
        trueStressStrainData(stressStrainIterator,2) = log(1 + inputFoamData(stressStrainIterator,2));
        trueStressStrainData(stressStrainIterator,3) = inputFoamData(stressStrainIterator,3) * (1 + inputFoamData(stressStrainIterator,2));
    end
end

%Correct the linear part of the curve to fit a linear elastic region
function [correctedLinearRegionData, foamDataStruct] = correctLinearRegion(inputFoamData, foamDataStruct)
    %foam data should be input as TIME, STRAIN, STRESS
    DEBUG = false;
    %Check Inputs
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 3), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME, STRAIN, STRESS");
    %Calculation
    %Rewrite the data betwene the strain offset and the yield point and correct for the strain offset
    %How many data points betwene yeild and zero
    numberOfDataPointsInLinearRegion = foamDataStruct.compressiveYieldStressIndex; 
    correctedLinearRegionData = zeros(inputDataSize(1),inputDataSize(2));

    %Rewrite the linear data
    strainStepSize = foamDataStruct.strainAtYield / (numberOfDataPointsInLinearRegion-1);

    for(rewriteModulusDataStreamIterator = 1:3)
        for(rewriteModulusDataIterator = 1:foamDataStruct.compressiveYieldStressIndex)
            switch(rewriteModulusDataStreamIterator)
                case 3
                    correctedLinearRegionData(rewriteModulusDataIterator,rewriteModulusDataStreamIterator) = foamDataStruct.modulus * strainStepSize * (rewriteModulusDataIterator-1);
                case 2
                    correctedLinearRegionData(rewriteModulusDataIterator,rewriteModulusDataStreamIterator) = strainStepSize * (rewriteModulusDataIterator-1);
                otherwise
                    correctedLinearRegionData(rewriteModulusDataIterator,rewriteModulusDataStreamIterator) = inputFoamData(rewriteModulusDataIterator,rewriteModulusDataStreamIterator);
            end
        end
    end
    
    %Copy over the data from the yield point onwards and correct for the strain offset
    for(rewriteModulusDataStreamIterator = 1:3) 
        for(rewriteModulusDataIterator = foamDataStruct.compressiveYieldStressIndex+1:length(inputFoamData))
            switch(rewriteModulusDataStreamIterator)
                case 2
                    correctedLinearRegionData(rewriteModulusDataIterator,rewriteModulusDataStreamIterator) = inputFoamData(rewriteModulusDataIterator,rewriteModulusDataStreamIterator)-foamDataStruct.strainOffsetValue;
                otherwise
                    correctedLinearRegionData(rewriteModulusDataIterator,rewriteModulusDataStreamIterator) = inputFoamData(rewriteModulusDataIterator,rewriteModulusDataStreamIterator);
            end
        end
    end
    %Do a calculation to work out the starting point then create the
    %offset points

    %Check the values are correct
    switch(DEBUG)
        case true
            %compressive yeild stress
            if(foamDataStruct.compressiveYieldStress == correctedLinearRegionData(foamDataStruct.compressiveYieldStressIndex,3))
                logDisp("The compressive yeild stress is correct",1,DEBUG);
            else
                logErr("The compressive yeild stress is wrong, check",-1,DEBUG);
            end
            %strain at yeild
            if((foamDataStruct.compressiveYieldStress / foamDataStruct.modulus) == correctedLinearRegionData(foamDataStruct.compressiveYieldStressIndex,2))
                logDisp("The compressive yeild strain is correct",1,DEBUG);
            else
                logErr("The compressive yeild strain is wrong, check",-1,DEBUG);
            end
        otherwise
    end

    %Rewrite indecies on the foam data struct becuse they have changed
    foamDataStruct.indexLocationUCS = foamDataStruct.indexLocationUCS;
    foamDataStruct.compressiveYieldStressIndex = foamDataStruct.compressiveYieldStressIndex;
    foamDataStruct.strainOffsetClosestIndex = 1;
    foamDataStruct.strainOffsetValue = 0;
    
end

%Cut out the plastic region and convert it to plastic strain for abaqus
function [outputPlasticStrainData] = convertStrainToAbaqusPlasticStrain(inputFoamData, foamDataStruct)
    %foam data should be input as TIME, STRAIN, STRESS
    %Check Inputs
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 3), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME, STRAIN, STRESS");
    %Calculation
    %Cut out the plastic region and convert it to plastic strain values
    outputPlasticStrainData = zeros(inputDataSize(1)-foamDataStruct.compressiveYieldStressIndex,inputDataSize(2));
    for(plasticRegionStreamSelector = 1:inputDataSize(2))
        for(plasticRegionDataSelector = foamDataStruct.compressiveYieldStressIndex:inputDataSize(1))
            switch(plasticRegionStreamSelector)
                case 2
                    %theoretically this is plastic strain et - el = et -
                    %(Stress/E)
                    outputPlasticStrainData(plasticRegionDataSelector-(foamDataStruct.compressiveYieldStressIndex-1),plasticRegionStreamSelector) = (inputFoamData(plasticRegionDataSelector,plasticRegionStreamSelector) - (inputFoamData(plasticRegionDataSelector,plasticRegionStreamSelector+1)/foamDataStruct.modulus));
                otherwise
                    outputPlasticStrainData(plasticRegionDataSelector-(foamDataStruct.compressiveYieldStressIndex-1),plasticRegionStreamSelector) = inputFoamData(plasticRegionDataSelector,plasticRegionStreamSelector);
            end
        end
    end
end

%Average the data accross the material 
function [stressAveragedData] = stressStrainCurveAverage(inputFoamData,strainStepSize)
    %Foam data should be input as TIME1,TIME2, STRAIN1, STRAIN2, STRESS1,
    %STRESS2
    %Equalize curve lengths
    inputFoamData = equaliseMaterialDataCurveLengths(inputFoamData);
    %Check Inputs
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && (mod(inputDataSize(2),3) == 0)), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME1, TIME2.., STRAIN1, STRAIN2.., STRESS1, STRESS2");
    %Calculation

    %Work out how many data streams there are and the average end strain
    inputNumberOfDataStreams = inputDataSize(2)/3;
    
    averageEndStrainData = 0;
    for(averageEndStrainIterator = (inputNumberOfDataStreams+1):(inputNumberOfDataStreams*2))
        averageEndStrainData = inputFoamData(end,averageEndStrainIterator)+averageEndStrainData;
    end
    averageEndStrainData = averageEndStrainData/inputNumberOfDataStreams;
    
    %Find the number of steps through the data for the strainStepSize
    assert((averageEndStrainData >= strainStepSize), "[STRAIN STEP SIZE ERROR]: THE STRAIN STEP SIZE IS LARGER THAN THE MAXIMUM STRAIN IN THE DATA");
    numberOfStepsInTheDataStream = ceil(averageEndStrainData/strainStepSize);
    
    %Average the data using "Square method" I.E.    1,3 = 1.5,3.5
    %                                               2,4
    assert((numberOfStepsInTheDataStream <= inputDataSize(1)), "[STRAIN STEP SIZE ERROR]: THE STEP SIZE IS SMALLER THAN STRAIN STEPS IN THE DATA");
    assert((numberOfStepsInTheDataStream < inputDataSize(1)), "[STRAIN STEP SIZE ERROR]: THE STEP SIZE IS TOO SMALL THERE WILL BE NO DATA AVERAGING");
    intermediateStressAveragedData = zeros(numberOfStepsInTheDataStream,3);
    numberOfDataPointsInEachStep = floor((inputDataSize(1))/numberOfStepsInTheDataStream);
    for(walkThroughDataAverageIterator = 1:numberOfStepsInTheDataStream)
        for(streamSelectorIterator = 1:3)
            tempAverage = 0;
            averagingCount = 0;
            for(numberOfDataStreamsIterator = 1:inputNumberOfDataStreams)
                for(numberOfDataPointsInEachStepIterator = 1:numberOfDataPointsInEachStep)
                    %((((walkThroughDataAverageIterator-1)*numberOfDataPointsInEachStep)+numberOfDataPointsInEachStepIterator))
                    switch(((((walkThroughDataAverageIterator-1)*numberOfDataPointsInEachStep)+numberOfDataPointsInEachStepIterator)==1)||((((walkThroughDataAverageIterator-1)*numberOfDataPointsInEachStep)+numberOfDataPointsInEachStepIterator)>=inputDataSize(1)))
                        case 1
                            switch(inputDataSize(1) <= (((walkThroughDataAverageIterator-1)*numberOfDataPointsInEachStep)+numberOfDataPointsInEachStepIterator))
                                case 1
                                    break;
                                otherwise
                                    continue;
                            end
                        otherwise
                            tempAverage = inputFoamData((((walkThroughDataAverageIterator-1)*numberOfDataPointsInEachStep)+numberOfDataPointsInEachStepIterator),(((streamSelectorIterator-1)*inputNumberOfDataStreams)+numberOfDataStreamsIterator)) + tempAverage;
                            averagingCount = averagingCount + 1;
                    end
                end
            end
            intermediateStressAveragedData(walkThroughDataAverageIterator,streamSelectorIterator) = (tempAverage / averagingCount);
        end
    end
    
    %First and last data points
    firstAndLastDataPoints = zeros(2,3);
    firstOrLastCount = 1;
    for(dataPointIterator = 1:(inputDataSize(1)-1):inputDataSize(1))
        for(streamSelector = 1:3)
            tempAverage = 0;
            averagingCount = 0;
            for(firstAndLastSummationIterator = 1:inputNumberOfDataStreams)
                tempAverage = inputFoamData((dataPointIterator),((streamSelector-1)*inputNumberOfDataStreams + firstAndLastSummationIterator)) + tempAverage;
                averagingCount = averagingCount + 1;
            end
            firstAndLastDataPoints(firstOrLastCount,streamSelector) = (tempAverage/averagingCount);
        end
        firstOrLastCount = firstOrLastCount + 1;
    end

    %Write all data to the same matrix
    stressAveragedData = zeros((length(intermediateStressAveragedData)+2),3);
    lineCount = 1;
    for(firstAndLastWriteIterator = 1:(length(stressAveragedData)-1):length(stressAveragedData))
        for(streamSelectorIterator = 1:3)
            stressAveragedData(firstAndLastWriteIterator,streamSelectorIterator) = firstAndLastDataPoints(lineCount,streamSelectorIterator);
        end
        lineCount = lineCount + 1;
    end
    for(writeDataIterator = 2:(length(intermediateStressAveragedData)+1))
        for(streamSelectorIterator = 1:3)
            stressAveragedData(writeDataIterator,streamSelectorIterator) = intermediateStressAveragedData(writeDataIterator-1,streamSelectorIterator);
        end
    end
end

%average material properties
function [averageDataVector] = averageMaterialProperties(inputDataVector)
    %foam data should be input as MODULUS, ULTIMATE COMPRESSIVE STRESS,
    %COMPRESSIVE YIELD STRESS
    %Check Inputs
    inputDataSize = size(inputDataVector);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 3), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: MODULUS, ULTIMATE COMPRESSIVE STRESS, COMPRESSIVE YIELD STRESS");
    %Calculation
    
    %Calculate average data
    averageDataVector = zeros(1,3);
    for(streamSelectorIterator = 1:inputDataSize(2))
        tempAverage = 0;
        averagingCount = 0;
        for(averagingDataIterator = 1:inputDataSize(1))
            tempAverage = inputDataVector(averagingDataIterator,streamSelectorIterator) + tempAverage;
            averagingCount = averagingCount + 1;
        end
        averageDataVector(1,streamSelectorIterator) = (tempAverage/averagingCount);
    end
end

%Find the location of the average yeild stress
function [outputMaterialStruct] = generateNewAverageMaterialDataStructOrriginalModulus(inputFoamData, inputAverageMaterialData)
    DEBUG = true;
    %Assumes input at TRUE STRESS STRAIN
    %foam data should be input as TIME, STRAIN, STRESS
    %Check Inputs
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 3), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME, STRAIN, STRESS");
    
    %foam data should be input as MODULUS, ULTIMATE COMPRESSIVE STRESS,
    %COMPRESSIVE YIELD STRESS
    %Check Inputs
    inputDataVectorSize = size(inputAverageMaterialData);
    sizeOfDataVectorSize = size(inputDataVectorSize);
    assert(((sizeOfDataVectorSize(2) == 2) && inputDataVectorSize(2) == 3), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: MODULUS, ULTIMATE COMPRESSIVE STRESS, COMPRESSIVE YIELD STRESS");
    %Calculation
    %New material data struct
    outputMaterialStruct = struct('modulus',0,'modulusIntercept',0,'ultimateCompressiveStress',0,'indexLocationUCS',0,'compressiveYieldStress',0,'compressiveYieldStressIndex',0);
    
    %find the location of the average compressive yeild stress
    yieldCompressiveStressIndex = 0;
    for(searchIterator = 1:inputDataSize(1))
        switch(inputFoamData(searchIterator,3) >= inputAverageMaterialData(1,3))
            case 1
                yieldCompressiveStressIndex = searchIterator;
                break;
            otherwise
                continue;
        end
    end
    outputMaterialStruct.compressiveYieldStressIndex = yieldCompressiveStressIndex;
    outputMaterialStruct.compressiveYieldStress = inputFoamData(yieldCompressiveStressIndex,3);
    
    %Look for UCS point
    ultimateCompressiveStressDataStruct = findUltimateCompressiveStressStruct(inputFoamData);
    outputMaterialStruct.ultimateCompressiveStress = ultimateCompressiveStressDataStruct.ultimateCompressiveStressValue;
    outputMaterialStruct.indexLocationUCS = ultimateCompressiveStressDataStruct.ultimateCompressiveStressValueLocation;
    
    %assign orriginal average modulus and offset
    outputMaterialStruct.modulus = inputAverageMaterialData(1,1);
    outputMaterialStruct.modulusIntercept = ((outputMaterialStruct.compressiveYieldStress - outputMaterialStruct.modulus * inputFoamData(yieldCompressiveStressIndex,2)));
   
    %Work out the strain axis intercept (Zero stress strain offset)
    strainOffsetValue = (-outputMaterialStruct.modulusIntercept)/outputMaterialStruct.modulus;
    
    %Find the strain index at the strain offset point
    for(strainOffsetSearchIterator = 1:length(inputFoamData))
        checkStrainComparitor = (inputFoamData(strainOffsetSearchIterator,2)>=strainOffsetValue);
        switch(checkStrainComparitor)
            case 1
                break;
            case 0
                continue;
        end
    end
    strainOffsetIndex = strainOffsetSearchIterator - 1;
    outputMaterialStruct.strainOffsetValue = strainOffsetValue;
    outputMaterialStruct.strainOffsetClosestIndex= strainOffsetIndex;

    logDisp("Checking the accuracy of the model",-1,DEBUG);
    if(outputMaterialStruct.ultimateCompressiveStress == inputFoamData(outputMaterialStruct.indexLocationUCS,3))
        logDisp("The ultimate compressive stress index is correct",1,DEBUG);
    else
        logErr("The ultimate compressive stress index is wrong, check",-1,DEBUG);
    end
end

%Generate a new output material struct
function [outputMaterialStruct] = generateNewAverageMaterialDataStructModifiedModulus(inputFoamData, inputAverageMaterialData)
    DEBUG = true;
    %Assumes input at TRUE STRESS STRAIN
    %foam data should be input as TIME, STRAIN, STRESS
    %Check Inputs
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 3), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME, STRAIN, STRESS");
    
    %foam data should be input as MODULUS, ULTIMATE COMPRESSIVE STRESS,
    %COMPRESSIVE YIELD STRESS
    %Check Inputs
    inputDataVectorSize = size(inputAverageMaterialData);
    sizeOfDataVecotorSize = size(inputDataVectorSize);
    assert(((sizeOfDataVecotorSize(2) == 2) && inputDataVectorSize(2) == 3), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: MODULUS, ULTIMATE COMPRESSIVE STRESS, COMPRESSIVE YIELD STRESS");
    %Calculation
    %New material data struct
    outputMaterialStruct = struct('modulus',0,'modulusIntercept',0,'ultimateCompressiveStress',0,'indexLocationUCS',0,'compressiveYieldStress',0,'compressiveYieldStressIndex',0);
    
    %find the location of the average compressive yeild stress
    yieldCompressiveStressIndex = 0;
    for(searchIterator = 1:inputDataSize(1))
        switch(inputFoamData(searchIterator,3) >= inputAverageMaterialData(1,3))
            case 1
                yieldCompressiveStressIndex = searchIterator;
                break;
            otherwise
                continue;
        end
    end
    outputMaterialStruct.compressiveYieldStressIndex = yieldCompressiveStressIndex;
    outputMaterialStruct.compressiveYieldStress = inputFoamData(yieldCompressiveStressIndex,3);
    
    %Look for UCS point
    ultimateCompressiveStressDataStruct = findUltimateCompressiveStressStruct(inputFoamData);
    outputMaterialStruct.ultimateCompressiveStress = ultimateCompressiveStressDataStruct.ultimateCompressiveStressValue;
    outputMaterialStruct.indexLocationUCS = ultimateCompressiveStressDataStruct.ultimateCompressiveStressValueLocation;
    
    %Correct the modulus data point using Linear Model Modulus
    modulusCalculationVector = zeros(outputMaterialStruct.compressiveYieldStressIndex,2);
    for(assignDataPointsIterator = 1:outputMaterialStruct.compressiveYieldStressIndex)
        for(streamSelectionIterator = 1:2)
            modulusCalculationVector(assignDataPointsIterator,streamSelectionIterator) = inputFoamData(assignDataPointsIterator,(streamSelectionIterator+1));
        end
    end
    fitTypeDesignation = 'linear';
    linearModel = fitlm(modulusCalculationVector(:,1),modulusCalculationVector(:,2),fitTypeDesignation);
    
    %Work out the strain axis intercept (Zero stress strain offset)
    strainOffsetValue = (-table2array(linearModel.Coefficients(1,1))/table2array(linearModel.Coefficients(2,1)));
    
    %Find the strain index at the strain offset point
    for(strainOffsetSearchIterator = 1:length(inputFoamData))
        checkStrainComparitor = (inputFoamData(strainOffsetSearchIterator,2)>=strainOffsetValue);
        switch(checkStrainComparitor)
            case 1
                break;
            case 0
                continue;
        end
    end
    strainOffsetIndex = strainOffsetSearchIterator - 1;
    
    %Write data to foam output data struct
    outputMaterialStruct.modulus = table2array(linearModel.Coefficients(2,1));
    outputMaterialStruct.modulusIntercept = table2array(linearModel.Coefficients(1,1));
    outputMaterialStruct.strainOffsetValue = strainOffsetValue;
    outputMaterialStruct.strainOffsetClosestIndex = strainOffsetIndex;

    logDisp("Checking the accuracy of the model",-1,DEBUG);
    if(outputMaterialStruct.ultimateCompressiveStress == inputFoamData(outputMaterialStruct.indexLocationUCS,3))
        logDisp("The ultimate compressive stress index is correct",1,DEBUG);
    else
        logErr("The ultimate compressive stress index is wrong, check",-1,DEBUG);
    end
end

%Cut the tail off the material data curve
function [lengthEqualisedOutputCurves] = equaliseMaterialDataCurveLengths(inputMaterialDataCurves)
    %Foam data should be input as TIME1,TIME2, STRAIN1, STRAIN2, STRESS1,
    %STRESS2
    
    %[ACTUALLY QUITE IMPORTANT]
    %Assumes the final consolidation part of the curves are the least
    %important so sacrifices them to presurve fidelity in the rest of the
    %curve
    
    %Check Inputs
    inputDataSize = size(inputMaterialDataCurves);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && (mod(inputDataSize(2),3) == 0)), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME1, TIME2.., STRAIN1, STRAIN2.., STRESS1, STRESS2");
    %Calculation
    %Work out how many data streams there are
    inputNumberOfDataStreams = inputDataSize(2)/3;
    assert(mod(inputNumberOfDataStreams,1)==0,"[DATA INPUT ERROR]: THE DATA INPUT IS NOT A MULTIPLE OF 3. PLEASE CHECK AND TRY AGAIN");
    
    %Work out which is the shortest data stream
    curveStartFlag = 0;
    streamCount = 0;
    for(searchLengthIterator = 2:inputDataSize(1))
        streamCountCheck = 0;
        for(streamSelector = 1:inputNumberOfDataStreams)
            switch(curveStartFlag == 0)
                case 1
                    switch(inputMaterialDataCurves(searchLengthIterator,streamSelector)==0)
                        case 1
                            continue;
                        case 0
                            curveStartFlag = 1;
                        otherwise
                            assert(1==0,"[CONTROL FLOW ERROR]: THE CURVE DATA POINT SHOULD BE ZERO OR NOT ZERO")
                    end
            end
            switch((inputMaterialDataCurves(searchLengthIterator,streamSelector)==0))
                case 1
                    break;
                case 0
                    switch(curveStartFlag==1)
                        case 1
                            streamCountCheck = streamCountCheck + 1;
                        case 0
                            continue;
                    end
            end
        end
        switch(streamCountCheck == inputNumberOfDataStreams)
            case 1
                streamCount = streamCount + 1;
            case 0
                break;
        end
    end
    
    %Now copy over all data other than the remaining portion of the
    lengthEqualisedOutputCurves = zeros(streamCount,inputDataSize(2));
    
    for(equalisedLengthStreamCopyIterator = 1:streamCount)
        for(streamCopyIterator = 1:inputDataSize(2))
            lengthEqualisedOutputCurves(equalisedLengthStreamCopyIterator,streamCopyIterator) = inputMaterialDataCurves(equalisedLengthStreamCopyIterator,streamCopyIterator);
        end
    end
end

%Produce equally distrbuted strain points, interpolate betwene data points
%linearly if nessesary.
function[outputDataPoints] = extractEquallyDistributedABAQUSPlasticStrainDataPoints(inputFoamData, inputFoamMaterialStruct, numberOfDataPoints)
    DEBUG = false;
    %Assumes input at TRUE STRESS STRAIN
    %foam data should be input as TIME, STRAIN, STRESS
    %Check Inputs
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 3), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME, STRAIN, STRESS");
    
    assert((inputDataSize(1) >= numberOfDataPoints),"[DATA LENGTH ERROR]: THE INPUT DATA CURVE HAS FEWER DATA POINTS THAN ARE REQUESTED IN OUTPUT");
    assert((3 <= numberOfDataPoints),"[DATA LENGTH ERROR]: THE REQUESTED NUMBER IF DATA POINTS IS INSUFFICIENT TO REPRESENT: YIELD, UCS, END");
    
    %Calculation
    arrayOfKeyLocations = zeros(2,3); %Store key locations beginning, ucs, and end: Stress and Strain
    
    arrayOfKeyLocations(1,1) = inputFoamData(1,3);%Beginning stress
    arrayOfKeyLocations(2,1) = 1;%Beginning index

    arrayOfKeyLocations(1,2) = inputFoamMaterialStruct.ultimateCompressiveStress;%Ultimate compressive stress
    arrayOfKeyLocations(2,2) = inputFoamMaterialStruct.indexLocationUCS; %Ultimate compressive strain

    %Find the maximum strain location
    maximumStrainLocation = findValueIndex(findMaximumValue(inputFoamData(:,2),1,10),inputFoamData(:,2),1);
    arrayOfKeyLocations(1,3) = inputFoamData(maximumStrainLocation,3);
    arrayOfKeyLocations(2,3) = maximumStrainLocation;


    %Check the values of UCS
    logDisp("Check the value of the UCS",-1 ,DEBUG);
    logDisp("Struct value",-1 ,DEBUG);
    logDisp(inputFoamMaterialStruct.ultimateCompressiveStress,-1 ,DEBUG);
    logDisp("Value at index location",-1 ,DEBUG);
    logDisp(inputFoamData(inputFoamMaterialStruct.indexLocationUCS,3),-1 ,DEBUG);

    %Find the strain range
    strainRange = inputFoamData(arrayOfKeyLocations(2,3),2) - inputFoamData(arrayOfKeyLocations(2,1),2);
    

    %Find the strain devisions
    strainRangeStep = strainRange / (numberOfDataPoints-1);

    outputDataPoints = zeros((numberOfDataPoints),2);
    
    %Create evenly spaced strain data
    for(strainDataAssignmentIterator = 1:size(outputDataPoints))
        outputDataPoints(strainDataAssignmentIterator,1) = strainRangeStep * (strainDataAssignmentIterator-1);
    end

    %Write out the data points and do linear interpolation betwene stress
    %points from strain location
    %check 0th special case and compensate for that
    outputDataPoints(1,2) = inputFoamData(1,3);

    %

    lastClosestStrainIndex = 1;
    for(outputDataPointIterator = 2:size(outputDataPoints))
        for(searchInputDataPointsIterator = lastClosestStrainIndex:size(inputFoamData,1))
            switch(inputFoamData(searchInputDataPointsIterator,2) >= outputDataPoints(outputDataPointIterator,1))
                case true
                    lastClosestStrainIndex = searchInputDataPointsIterator;
                    outputDataPoints(outputDataPointIterator,2) = linearInterpolation(inputFoamData((searchInputDataPointsIterator-1):searchInputDataPointsIterator,2:3),outputDataPoints(outputDataPointIterator,1));
                    break;
                otherwise
                    continue;
            end
        end
    end

    fprintf("\nData points desired: %d\n",numberOfDataPoints);
    fprintf("Data points in data: %d\n",length(outputDataPoints));
    fprintf("\n");
end

%[ERROR]
%Select an equally distributed 20 points from the material data point 
% (10 before UCS, 6 in 0.75 msx strain, 4 last 0.25 max strain)
%THIS FUNCTION WILL CAUSE AN ERROR IN ABAQUS AT THIS TIME
%ABAQUS REQUIRES EQUALLY DISTRIBUTED DATA POINTS
%[ERROR]
function [outputDataPoints] = extractPrioritizedABAQUSPlasticStrainDataPoints(inputFoamData ,numberOfDataPoints)
    %Assumes input at TRUE STRESS STRAIN
    %foam data should be input as TIME, STRAIN, STRESS
    %Check Inputs
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 3), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME, STRAIN, STRESS");
    
    assert((inputDataSize(1) >= numberOfDataPoints),"[DATA LENGTH ERROR]: THE INPUT DATA CURVE HAS FEWER DATA POINTS THAN ARE REQUESTED IN OUTPUT");
    assert((3 <= numberOfDataPoints),"[DATA LENGTH ERROR]: THE REQUESTED NUMBER IF DATA POINTS IS INSUFFICIENT TO REPRESENT: YIELD, UCS, END");
    %Calculation

    arrayOfLocationRatios = zeros(3,1);%from start to ucs, 1; from ucs to 75% strain, 2; from 75% strain
    arrayOfLocationRatios(1,1) = 0.5;
    arrayOfLocationRatios(2,1) = 0.3;
    arrayOfLocationRatios(3,1) = 1-(arrayOfLocationRatios(1,1) + arrayOfLocationRatios(2,1));

    numberMinusKeyPoints = numberOfDataPoints - 3;
    arrayOfDataPointsBetweenKeyLocations = zeros(3,1);
    arrayOfDataPointsBetweenKeyLocations(1,1) = ceil(arrayOfLocationRatios(1) * numberMinusKeyPoints); 
    arrayOfDataPointsBetweenKeyLocations(3,1) = ceil(arrayOfLocationRatios(3) * numberMinusKeyPoints); 
    arrayOfDataPointsBetweenKeyLocations(2,1) = numberMinusKeyPoints - (arrayOfDataPointsBetweenKeyLocations(3,1) + arrayOfDataPointsBetweenKeyLocations(1,1));

    idealDistributionOfDataPoints = arrayOfDataPointsBetweenKeyLocations;
    %Find index of key points
    arrayOfKeyIndexPoints = zeros(4,1); %{YeildPoint,1; UCS,2; 75%strain,3; end,4}
    arrayOfKeyIndexPoints(1,1) = 1;

    %Search for the largest strain point and assign this as the end strain
    %point
    endOfStrainDataFlag = 0;
    maximumStrainIndex = 0;
    foundMaximumStrainFlag = 0;
    numberOfFutureDataPointsToSearch = 50;
    comparisonVariable = findMaximumValue(inputFoamData(:,2),1,numberOfFutureDataPointsToSearch);
    comparisonVariableIndex = findValueIndex(comparisonVariable,inputFoamData(:,2),1,"all");

%   FEELTING HAVE LEFT OUT THE DATA COMBING FUNCTION IN THIS BUT IT SHOULD
%   BE IN A SEPERATE FUNCTION ANYWAY

    arrayOfKeyIndexPoints(4,1) = findValueIndex(findMaximumValue(inputFoamData(:,2),1,numberOfFutureDataPointsToSearch),inputFoamData(:,2),1);
    assert(((arrayOfKeyIndexPoints(4,1)-3) >= numberOfDataPoints),"[DATA LENGTH ERROR]: THE RECALCULATED LENGTH STRAIN CURVE HAS FEWER DATA POINTS THAN ARE REQUESTED IN OUTPUT");

    ultimateCompressionStress = inputFoamData(1,3);
    ultimateCompressiveStressDataStruct = findUltimateCompressiveStressStruct(inputFoamData);
    arrayOfKeyIndexPoints(2,1) = ultimateCompressiveStressDataStruct.ultimateCompressiveStressValueLocation;

    %Find index location for 75% strain location
    strain75OfMaxValue = inputFoamData(end,2) * 0.75;
    for(searchFor75StrainLocationIterator = arrayOfKeyIndexPoints(2,1):length(inputFoamData))
        switch(inputFoamData(searchFor75StrainLocationIterator,2) >= strain75OfMaxValue)
            case 1 
                arrayOfKeyIndexPoints(3,1) = searchFor75StrainLocationIterator;
                break;
            otherwise
                continue;
        end
    end

    %If the data points do not exist - redistribute to the other sections
    %based on their relative ratio
    %Then check the other sections to see if data points exist and so of
    for(numberOfSectionsIterator = 1:length(arrayOfKeyIndexPoints)-1)
        switch((arrayOfKeyIndexPoints(numberOfSectionsIterator+1,1) - arrayOfKeyIndexPoints(numberOfSectionsIterator,1)-1)>=arrayOfDataPointsBetweenKeyLocations(numberOfSectionsIterator,1))
            case 0
                arrayOfDataPointsBetweenKeyLocations(numberOfSectionsIterator+1,1) = (arrayOfDataPointsBetweenKeyLocations(numberOfSectionsIterator,1) - (arrayOfKeyIndexPoints(numberOfSectionsIterator+1,1) - arrayOfKeyIndexPoints(numberOfSectionsIterator,1)-1)) + arrayOfDataPointsBetweenKeyLocations(numberOfSectionsIterator+1,1);
                arrayOfDataPointsBetweenKeyLocations(numberOfSectionsIterator, 1) = (arrayOfKeyIndexPoints(numberOfSectionsIterator+1,1) - arrayOfKeyIndexPoints(numberOfSectionsIterator,1)-1);
            case 1
                continue;
        end
    end

    %Work out how to take steps through the data
    arrayOfStepsBetweenDataPointsInSections = zeros(length(arrayOfDataPointsBetweenKeyLocations),2);
    for(stepCalculateIterator = 1 : length(arrayOfStepsBetweenDataPointsInSections))
        arrayOfStepsBetweenDataPointsInSections(stepCalculateIterator,1) = ceil((arrayOfKeyIndexPoints(stepCalculateIterator+1,1) - arrayOfKeyIndexPoints(stepCalculateIterator,1))/(arrayOfDataPointsBetweenKeyLocations(stepCalculateIterator,1)+1));
        switch((1 > mod((arrayOfKeyIndexPoints(stepCalculateIterator+1,1) - arrayOfKeyIndexPoints(stepCalculateIterator,1))/(arrayOfDataPointsBetweenKeyLocations(stepCalculateIterator,1)+1),1)) && (0 < mod((arrayOfKeyIndexPoints(stepCalculateIterator+1,1) - arrayOfKeyIndexPoints(stepCalculateIterator,1))/(arrayOfDataPointsBetweenKeyLocations(stepCalculateIterator,1)+1),1)))
            case 1
                arrayOfStepsBetweenDataPointsInSections(stepCalculateIterator,2) = ceil(1/mod((arrayOfKeyIndexPoints(stepCalculateIterator+1,1) - arrayOfKeyIndexPoints(stepCalculateIterator,1))/(arrayOfDataPointsBetweenKeyLocations(stepCalculateIterator,1)+1),1));
            otherwise
                continue;
        end
    end

    %Assign the material indecies to an arrray to allow easy assignment to
    %the output array
    arrayOfIndecies = zeros(numberOfDataPoints,1);
    assignmentCounter = 1;
    for(sectionIterator = 1:length(arrayOfDataPointsBetweenKeyLocations))
        extraDataPointInsertCounter = 0;
        for(itterateThroughKeyLocationGapIterator = arrayOfKeyIndexPoints(sectionIterator,1):arrayOfStepsBetweenDataPointsInSections(sectionIterator,1):arrayOfKeyIndexPoints(sectionIterator+1,1)-1)
            arrayOfIndecies(assignmentCounter,1) = itterateThroughKeyLocationGapIterator;
            switch((arrayOfStepsBetweenDataPointsInSections(sectionIterator,2)~=0)&&(arrayOfStepsBetweenDataPointsInSections(sectionIterator,2)==extraDataPointInsertCounter))
                case 1
                    assignmentCounter = assignmentCounter + 1;
                    arrayOfIndecies(assignmentCounter,1) = itterateThroughKeyLocationGapIterator + ceil(arrayOfStepsBetweenDataPointsInSections(sectionIterator,1)/2);
                    extraDataPointInsertCounter = 0;
                case 0
                    extraDataPointInsertCounter = extraDataPointInsertCounter + 1;
            end
            assignmentCounter = assignmentCounter + 1;
        end
    end

    %Assign last data point
    arrayOfIndecies(end,1) = arrayOfKeyIndexPoints(end,1);

    %Search for zeros in the index data
    zeroCount = 0;
    for(searchForZerosInIndexList = 1:length(arrayOfIndecies))
        switch(arrayOfIndecies(searchForZerosInIndexList,1) == 0)
            case 1
                zeroCount = zeroCount + 1;
            otherwise
                continue;
        end
    end
    checkArrayOfIndecies = zeros(length(arrayOfIndecies) - zeroCount,1);
    nonZeroCount = 0;
    for(assignNonZerosInIndexList = 1:length(arrayOfIndecies))
        switch(arrayOfIndecies(assignNonZerosInIndexList,1) ~= 0)
            case 1
                nonZeroCount = nonZeroCount + 1;
                checkArrayOfIndecies(nonZeroCount,1) = arrayOfIndecies(assignNonZerosInIndexList,1);
            otherwise
                continue;
        end
    end

    %Assign data points to the output array
    outputDataPoints = zeros(length(checkArrayOfIndecies),2);
    for(assignmentIterator = 1:length(outputDataPoints))
        for(streamSelectionIterator = 1:2)
            outputDataPoints(assignmentIterator,streamSelectionIterator) = inputFoamData(checkArrayOfIndecies(assignmentIterator,1),streamSelectionIterator+1);
        end
    end
    outputDataPoints(1,1) = 0;

    fprintf("\nData points desired: %d\n",numberOfDataPoints);
    fprintf("Data points in data: %d\n",length(outputDataPoints));
    fprintf("\n");
    for(printIterator = 1:length(arrayOfLocationRatios))
        fprintf("\nDesired Data split %d: %0.3f",printIterator,(idealDistributionOfDataPoints(printIterator,1)/numberOfDataPoints));
    end
    fprintf("\n");fprintf("\n");
    for(printIterator = 1:length(arrayOfLocationRatios))
        fprintf("\nActual Data split %d: %0.3f",printIterator,(arrayOfDataPointsBetweenKeyLocations(printIterator,1)/length(outputDataPoints)));
    end
    fprintf("\n");fprintf("\n");
end

%Fracture energy function for plastic region
function [outputMaterialDataStruct] = getABAQUSPlasticFractureEnergy(inputFoamData, inputMateralDataStruct)
    %Assumes input at TRUE STRESS STRAIN
    %foam data should be input as STRESS, STRAIN
    %Check Inputs
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 2), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME, STRAIN, STRESS");
    
    %function trapz(vel)
    outputMaterialDataStruct = inputMateralDataStruct;
    outputMaterialDataStruct.plasticToughnessEnergy = trapz(inputFoamData(:,1),inputFoamData(:,2));
end

%Find the index or array of indecies that corrisponds to a desired value
function [outputIndex] = findValueIndex(varargin)
    DEBUG = false;
    %Input data is assumed to be in format:

    %SEARCH_FOR_VALUE, INPUT_FOAM_DATA, STREAM_TO_SEARCH_ALONG - (assumes search whole range)
    
    %Assign values and check inputs are alligend with the function expected
    %inputs
    assert(((nargin == 3)||(nargin == 4)),"[NUMBER OF INPUTS ERROR]: PLEASE INPUT DATA IN EITHER: SEARCH_FOR_VALUE, INPUT_FOAM_DATA, STREAM_TO_SEARCH_ALONG - (assumes search range for single value) OR SEARCH_FOR_VALUE, INPUT_FOAM_DATA, STREAM_TO_SEARCH_ALONG, MULTIPLE_DESIGNATOR - (finds the index location of all values equal to search)");
    
    %Check types
    searchForValue = varargin{1,1};
    assert(isa(searchForValue,'double'),"[INPUT TYPE ERROR]: THE INPUT VALUE TO SEARCH FOR MUST BE TYPE DOUBLE");
    inputFoamData = varargin{1,2};
    assert(isa(inputFoamData,'double'),"[INPUT TYPE ERROR]: THE INPUT MATERIAL DATA MUST BE TYPE DOUBLE");
    searchStreamIndex = varargin{1,3};
    assert(isa(searchStreamIndex,'double'),"[INPUT TYPE ERROR]: THE INPUT INDEX OF THE STREAM TO SEARCH IN MUST BE TYPE DOUBLE");
    
    multipleValuesFlag = false;

    switch(nargin)
        case 3
            logDisp("Three inputs",-1,DEBUG);
            multipleValuesFlag = false;
        otherwise
            logDisp("Four inputs",-1,DEBUG);
            multipleSearchValueString = varargin{1,4};
            assert(isa(multipleSearchValueString,'string'),"[INPUT TYPE ERROR]: THE MULTIPLE INDEX SEARCH DESIGNATOR SHOULD BE OF TYPE STRING");
            assert((multipleSearchValueString == "all"), "[INPUT ARGUMENT ERROR]: THE MULTIPLE INDEX SEARCH DESIGNATOR SHOULD BE: all");
            multipleValuesFlag = true;
    end

    %Check that the data is referenceable
    %The input search foam material data is assumed to be TIME, STRAIN, STRESS
    %This informs how the search will be conducted
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(inputDataSize(2)>=searchStreamIndex,"[STREAM SELECT ERROR]: THE SEARCH STREAM INDEX IS OUTSIDE OF THE INPUT FOAM INDEX BOUNDS");

    logDisp("Checks complete",-1,DEBUG);

    foundValueFlag = false;
    valueIndeciesVector  = [];
    indexCounter = 1;

    for(searchIterator = 1:inputDataSize(1))
        logDisp("Searching index",searchIterator,DEBUG);
        switch(searchForValue == inputFoamData(searchIterator,searchStreamIndex))
            case true
                logWarn("Found corrisiponding index",searchIterator,DEBUG);
                foundValueFlag = true; %Could come back to this and allow multiple value searches
                valueIndeciesVector(1,indexCounter) = searchIterator;
                indexCounter = indexCounter + 1;
            otherwise
                continue;
        end
        %Check there are not multiple values
        switch((foundValueFlag == true) && (multipleValuesFlag == false))
            case true
                break;
            otherwise
                continue;
        end
    end
    if(foundValueFlag == false)
        outputIndex = -1;
    elseif(multipleValuesFlag == false)
        outputIndex = valueIndeciesVector(1,1);
    else
        outputIndex = valueIndeciesVector;
    end
end

%Find the maximum value in an array and check future in the data for a higher value
function [outputValue] = findMaximumValue(varargin)
    DEBUG = false;
    %Will find the first largest number unless given a future search length
    %If given a future search length will find the first largest value in
    %including the future length
    %if equally large values are found will out put the first one
    %Input data is assumed to be in format:

    %INPUT_FOAM_DATA, STREAM_TO_SEARCH_ALONG, FUTURE_SEARCH_LENGTH 
    
    %Assign values and check inputs are alligend with the function expected
    %inputs
    assert(((nargin == 2)||(nargin == 3)),"[NUMBER OF INPUTS ERROR]: PLEASE INPUT DATA IN: INPUT_FOAM_DATA, STREAM_TO_SEARCH_ALONG OR INPUT_FOAM_DATA, STREAM_TO_SEARCH_ALONG, FUTURE_SEARCH_LENGTH");
    
    %Check types
    inputFoamData = varargin{1,1};
    assert(isa(inputFoamData,'double'),"[INPUT TYPE ERROR]: THE INPUT MATERIAL DATA MUST BE TYPE DOUBLE");
    searchStreamIndex = varargin{1,2};
    assert(isa(searchStreamIndex,'double'),"[INPUT TYPE ERROR]: THE INPUT INDEX OF THE STREAM TO SEARCH IN MUST BE TYPE DOUBLE");

    switch(nargin)
        case 2
            futureSearchLength = 0;
        otherwise
            futureSearchLength = varargin{1,3};
    end

    assert(isa(futureSearchLength,'double'),"[INPUT TYPE ERROR]: THE INPUT SEARCH FUTURE LENGTH MUST BE OF TYPE DOUBLE");

    %Check that the data is referenceable
    %The input search foam material data is assumed to be TIME, STRAIN, STRESS
    %This informs how the search will be conducted
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((inputDataSize(2)>=searchStreamIndex)&&(0<searchStreamIndex)),"[STREAM SELECT ERROR]: THE SEARCH STREAM INDEX IS OUTSIDE OF THE INPUT FOAM INDEX BOUNDS");

    assert(futureSearchLength < inputDataSize(1),"[INPUT ERROR]: FUTURE SEARCH MUST BE SHORTER THAN DATA STREAM");

    logDisp("The index seems to be ok",-1,DEBUG);

    foundMaximumValueFlag = false;
    outputIndex = 0;
    
    switch(inputDataSize(1)>1)
        case true
            logDisp("Going Into Search",-1,DEBUG);
            for(searchIterator = 1:(inputDataSize(1)-futureSearchLength-1))
                logDisp("Search Index",searchIterator,DEBUG);
                switch((inputFoamData(searchIterator,searchStreamIndex) >= inputFoamData(searchIterator+1,searchStreamIndex)))
                    case true
                        foundMaximumValueFlag = true;
                        logDisp("Going Into Future Search",-1,DEBUG);
                        for(futureSearchIterator = (searchIterator+1):(searchIterator+futureSearchLength))
                            usedFutureSearch = true;
                            logDisp("Future Search Index",futureSearchIterator,DEBUG);
                            switch(inputFoamData(searchIterator,searchStreamIndex) < inputFoamData(futureSearchIterator,searchStreamIndex))
                                case true
                                    foundMaximumValueFlag = false;
                                    break;
                                otherwise
                                    foundMaximumValueFlag = true;
                                    continue;
                            end
                        end
                        switch(foundMaximumValueFlag == true)
                            case true
                                outputIndex = searchIterator;
                            otherwise
                                continue;
                        end
                    otherwise
                        continue;
                end
                
                switch(foundMaximumValueFlag == true)
                    case true
                        break;
                    otherwise
                        continue;
                end
            end
        
            %Check last and second to last value and that the first value
            %is not the largest
            switch((inputFoamData(1,searchStreamIndex) >= inputFoamData(end-1,searchStreamIndex)) && (inputFoamData(1,searchStreamIndex) >= inputFoamData(end-2,searchStreamIndex)))
                case true
                    foundMaximumValueFlag = true;
                    outputIndex = 1;
                    logDisp("First value selected", -1, DEBUG);
                case false
                    switch((foundMaximumValueFlag == false)  && (inputFoamData(end,searchStreamIndex) >= inputFoamData(end-1,searchStreamIndex)))
                        case true
                            foundMaximumValueFlag = true;
                            outputIndex = inputDataSize(1);
                            logDisp("Last value selected", -1, DEBUG);
                        case false
                            foundMaximumValueFlag = true;
                            outputIndex = inputDataSize(1) - 1;
                            logDisp("Second to last value selected", -1, DEBUG);
                    end
            end

        otherwise
            foundMaximumValueFlag = true;
            outputIndex = 1;
            logDisp("Only one value input",1,DEBUG);
    end
    
    %Guard check on value to make sure the algorithm ran
    assert((foundMaximumValueFlag == true),"[ALGORITHM ERROR]: UNABLE TO FIND A MAXIMUM VALUE WITH THE PARAMETERS ENTERED, PLEASE ALTER AND TRY AGIAN");
    outputValue = inputFoamData(outputIndex,searchStreamIndex);
end

%linear interpolation function
function predictedY3Coordinate = linearInterpolation(varargin)
    %Expected input [[x1,y1];[x2,y2]],x3
    assert(nargin == 2, "[INPUT ERROR]: PLEASE INPUT [[X1,Y1];[X2,Y2]], X3; OUTPUT WILL BE Y3");
    
    vectorOfInputs = zeros(2,2);
    vectorOfInputs = varargin{1,1};
    strainAtDesiredPoint = varargin{1,2};

    %Calculation
    gradientOfLine = (vectorOfInputs(1,2) - vectorOfInputs(2,2)) / (vectorOfInputs(1,1) - vectorOfInputs(2,1));
    yInterceptOfLine = vectorOfInputs(1,2) - (vectorOfInputs(1,1) * gradientOfLine);

    %Output y = mx + c value
    predictedY3Coordinate = (gradientOfLine * strainAtDesiredPoint) + yInterceptOfLine;
end

%Feature detection functions%
%Detection of dwell feature of the Dartec
function [dartec_dwell_point,number_of_steps,time_step] = dartec_dwell_detect_function(input_data, time_data_stream_select, displacement_data_stream_select, test_speed)
    time_step = (((input_data(10,time_data_stream_select)-(input_data(1,time_data_stream_select))))/10);
    number_of_steps = (0.011/test_speed)/time_step; %value derived from round to upper .001 from maximum variation betwene two points over 15138 points in a dwell phase
    number_of_steps = ceil(number_of_steps)*10;
    distance_should_travel = (time_step * number_of_steps) * test_speed;
    for(itteration_variable = 1:(length(input_data)-number_of_steps))
        comparitor_variable = input_data((itteration_variable+number_of_steps),displacement_data_stream_select)-input_data(itteration_variable,displacement_data_stream_select);
        if((comparitor_variable <= (distance_should_travel))&&(1))
            dartec_dwell_point = itteration_variable;
            break;
        end
    end
end

%Value detection function
function value_location = value_detection_function(input_data, data_stream_select,find_value)
    value_location = 0;
    for(itterator_value = 1:length(input_data))
        if(find_value == input_data(itterator_value,data_stream_select))
            value_location = itterator_value;
            break;
        end
    end
    if((itterator_value == length(input_data)) && (value_location == 0))
        fprintf("The %.3f value was not found",find_value);
    end
end

%Find the material UCS point (Defined as either the local maximum or the
% local point where the second derivative changes direction
function [outputDataStruct] = findUltimateCompressiveStressStruct(inputFoamData)
    TIME = false;
    switch(TIME)
        case 1
            tic;
    end
    %DebugFunction
    DEBUG = false;
    %foam data should be input as TIME, STRAIN, STRESS
    %Check Inputs
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 3), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME, STRAIN, STRESS");
    %Calculation
    

    %Thought process
    %Write data to be fitted to a holding vector (Prealocated for speed)
    %This will half step with the fitting step
    %Look over the whole data stream find the "Section" of the data where the
    %asymptote happens

    %Walk along the curve and find the inflection point specifically the
    %point where the second derivative goes from positive to negative
    %Find maximum stress at that location

    %Initialise the search
    %These parameters have been tested on a number of different densities
    %and work
    numberOfSectionsToSplitTheDataInto = 20;
    minimumLengthOfDataSearch = 5;
    vectorOfCurveFits = zeros(numberOfSectionsToSplitTheDataInto,6);
    searchAreaStart = 1;
    searchAreaEnd = inputDataSize(1);
    searchLength = findSearchLength(searchAreaStart, searchAreaEnd, numberOfSectionsToSplitTheDataInto);
    switch(DEBUG)
        case 1
            n = 1;
            runCount = 1;
    end
    while(searchLength >= minimumLengthOfDataSearch)
        startSearch = searchAreaStart;
        endSearch = startSearch+searchLength;
        for(numberOfTimesToRunIterator = 1:numberOfSectionsToSplitTheDataInto)
            %Write fitting data to a test vector
            outOfSequenceIterator = 1;
            checkFittingVector = zeros(searchLength,3);
            for(writingLengthIterator = startSearch:endSearch)
                for(writingWidthPoints = 1:3)
                    checkFittingVector(outOfSequenceIterator,writingWidthPoints) = inputFoamData(writingLengthIterator,writingWidthPoints);
                end
                outOfSequenceIterator = outOfSequenceIterator+1;
            end
            %Remap Start and end points
            vectorOfCurveFits(numberOfTimesToRunIterator,4) = startSearch;
            vectorOfCurveFits(numberOfTimesToRunIterator,5) = endSearch;
            startSearch = endSearch;
            endSearch = startSearch + searchLength;

            %Fitting functions to the data
            %Linear Fit
            fitType = fittype("poly1");
            fitData = fit(checkFittingVector(:,2),checkFittingVector(:,3),fitType);
            switch(DEBUG)
                case 1
                    figure(n)
                    plot(fitData,checkFittingVector(:,2),checkFittingVector(:,3));
                    hold on
                    title(strcat("RunCount ",num2str(runCount)))
            end
            
    
            %Gather vector of fits
            vectorOfCurveFits(numberOfTimesToRunIterator,1) = fitData.p1;
            vectorOfCurveFits(numberOfTimesToRunIterator,2) = fitData.p2;
            vectorOfCurveFits(numberOfTimesToRunIterator,3) = 0;


        end
        switch(DEBUG)
            case 1
                n = n + 1;
        end

        %Find Gradient betwene points in the vectorOfCurveFits
        for(gradientCalculationIterator = 1:(size(vectorOfCurveFits,1)-1))
            vectorOfCurveFits(gradientCalculationIterator,6) = (vectorOfCurveFits(gradientCalculationIterator+1,1)-vectorOfCurveFits(gradientCalculationIterator,1));
        end
        %Find the angle betwene gradient 

        %Look for a greater then more neagtive then greater again take the
        %middle of that pattern
        for(searchForSignChangeIterator = 1:(size(vectorOfCurveFits,1)-1))
            if((vectorOfCurveFits(searchForSignChangeIterator,6) < 0) && (vectorOfCurveFits(searchForSignChangeIterator+1,6) > 0))
                if(searchForSignChangeIterator == 1)
                else
                    newSearchRegionStart = searchForSignChangeIterator-1;
                    newSearchRegionEnd = searchForSignChangeIterator+1;
                    break;
                end
            end
        end
        searchAreaStart = vectorOfCurveFits(newSearchRegionStart,4);
        searchAreaEnd = vectorOfCurveFits(newSearchRegionEnd,5);
        startSearch = searchAreaStart;

        switch(DEBUG)
            case 1
                figure(n)
                n = n + 1;
                plot(vectorOfCurveFits(:,1));
                currentGradientLocation = newSearchRegionStart;
                title(strcat("RunCount ",num2str(runCount)," - ChoseLocation ",num2str(currentGradientLocation)));
                runCount = runCount + 1;
        end
        searchLength = findSearchLength(searchAreaStart, searchAreaEnd, numberOfSectionsToSplitTheDataInto);
    end

    %Then look for the maximum Value in the range
    stressStreamIndex = 3;
    futureSearchSize = ceil((1/100)*(searchAreaEnd - searchAreaStart));
    outputDataStruct.ultimateCompressiveStressValue = findMaximumValue(inputFoamData(searchAreaStart:searchAreaEnd,stressStreamIndex),1,futureSearchSize);
    outputDataStruct.ultimateCompressiveStressValueLocation = findValueIndex(outputDataStruct.ultimateCompressiveStressValue,inputFoamData(:,stressStreamIndex),1);

    logDisp("Checking the accuracy of the model",-1,DEBUG);
    if(outputDataStruct.ultimateCompressiveStressValue == inputFoamData(outputDataStruct.ultimateCompressiveStressValueLocation,3))
        logDisp("The ultimate compressive stress index is correct",1,DEBUG);
    else
        logErr("The ultimate compressive stress index is wrong, check",-1,DEBUG);
    end

    %More debug features
    switch(DEBUG)
        case 1
            if(runCount > 1)
                endString = "S";
            else
                endString = " ";
            end
            fprintf("THE UCS SEARCHING ALGORTHYM RAN %d TIME%s\n",runCount,endString);
    end
    switch(TIME)
        case 1
            toc;
    end
end

function [outputSearchLength] = findSearchLength(inputSearchStart, inputSearchEnd, inputNumberOfSectionsToSplitTheDataInto)
    outputSearchLength = floor((inputSearchEnd - inputSearchStart)/inputNumberOfSectionsToSplitTheDataInto);
end

%Data output functions
%Write generic foam data to file |||||||||||||||||||||||||||||||||||||Need
%to fix at some points
function [] = writeToFileFoamData(orriginalFileName,inputData,timeLocation,displacementLocation,loadLocation)
    fileName = strcat("CutFoamData",orriginalFileName,".txt");
    fileID = fopen(fileName,"w");
    switch(fileID)
        case(-1)
            fprintf(strcat("File: ",fileName," ", "opening failed\n"));
            return;
        case(fileID)
            fprintf(strcat("File: ",fileName," ", "opened successfully\n"));
        otherwise
            assert(1==0,"[FILE OPEN ID SWITCH ERROR]: OUTSIDE CASE STRUCTURE");
    end
    
    fprintf(fileID,"TimeData\tStrainData\tStressData\n");
    for(itterator = 1:length(inputData))
        fprintf(fileID,"%.6f\t%.6f\t%.6f\n",inputData(itterator,timeLocation),inputData(itterator,displacementLocation),inputData(itterator,loadLocation));
    end

    file_close_comparator = fclose(fileID);
    switch(file_close_comparator)
        case(-1)
            fprintf(strcat("File: ",fileName," ", "closing failed\n"));
            return;
        case(0)
            fprintf(strcat("File: ",fileName," ", "closed successfully\n"));
        otherwise
            assert(1==0,"[FILE CLOSE ID SWITCH ERROR]: OUTSIDE CASE STRUCTURE");
    end
end

%Write ABAQUS input data to text file for transport to other people
function [] = writeABAQUSDataToFile(finalFileName,inputFoamData, inputMaterialDataStruct)
    %Assumes input at TRUE STRESS STRAIN
    %foam data should be input as STRESS, STRAIN
    %Check Inputs
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 2), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME, STRAIN, STRESS");
    
    %function
    fileName = strcat(finalFileName,".dat");
    fileID = fopen(fileName,"w");
    switch(fileID)
        case(-1)
            fprintf(strcat("File: ",fileName," ", "opening failed\n"));
            return;
        case(fileID)
            fprintf(strcat("File: ",fileName," ", "opened successfully\n"));
        otherwise
            assert(1==0,"[FILE OPEN ID SWITCH ERROR]: OUTSIDE CASE STRUCTURE");
    end
    
    fprintf(fileID,"BASIC MATERIAL DATA\n\n");
    fprintf(fileID,"MODULUS (Pa): \t\t\t\t%0.6f\n",inputMaterialDataStruct.modulus);
    fprintf(fileID,"COMPRESSIVE YIELD STRESS (Pa): \t\t%0.6f\n",inputMaterialDataStruct.compressiveYieldStress);
    fprintf(fileID,"ULTIMATE COMPRESSIVE STRESS (Pa): \t%0.6f\n",inputMaterialDataStruct.ultimateCompressiveStress);
    fprintf(fileID,"MATERIAL MAXIMUM STRAIN: \t\t%0.6f\n",inputFoamData(end,1));
    fprintf(fileID,"MATERIAL PLASTIC TOUGHNESS ENERGY: \t%0.6f\n",inputMaterialDataStruct.plasticToughnessEnergy);

    fprintf(fileID,"\n\nPLASTIC STRESS STRAIN DATA\n\n");
    fprintf(fileID,"StressData (Pa)\tStrainData\n\n");
    for(fileWriteIterator = 1:inputDataSize(1))
        for(streamIterator = inputDataSize(2):-1:1)
            fprintf(fileID,"%.6f\t",inputFoamData(fileWriteIterator,streamIterator));
        end
        fprintf(fileID,"\n");
    end

    file_close_comparator = fclose(fileID);
    switch(file_close_comparator)
        case(-1)
            fprintf(strcat("File: ",fileName," ", "closing failed\n"));
            return;
        case(0)
            fprintf(strcat("File: ",fileName," ", "closed successfully\n"));
        otherwise
            assert(1==0,"[FILE CLOSE ID SWITCH ERROR]: OUTSIDE CASE STRUCTURE");
    end
end

%Write ABAQUS input data to text file in python format
function [] = writeABAQUSPythonDataToFile(finalFileName,inputFoamData, inputMaterialDataStruct)
    %Assumes input at TRUE STRESS STRAIN
    %foam data should be input as STRESS, STRAIN
    %Check Inputs
    inputDataSize = size(inputFoamData);
    sizeOfDataSize = size(inputDataSize);
    assert(((sizeOfDataSize(2) == 2) && inputDataSize(2) == 2), "[DATA FORMAT ERROR]: THE DATA IS IN THE WRONG FORMAT. PLEASE CHECK INPUT DATA. FORMAT SHOULD BE: TIME, STRAIN, STRESS");
    
    %function
    fileName = strcat("PythonInput",finalFileName,".dat");
    fileID = fopen(fileName,"w");
    switch(fileID)
        case(-1)
            fprintf(strcat("File: ",fileName," ", "opening failed\n"));
            return;
        case(fileID)
            fprintf(strcat("File: ",fileName," ", "opened successfully\n"));
        otherwise
            assert(1==0,"[FILE OPEN ID SWITCH ERROR]: OUTSIDE CASE STRUCTURE");
    end
    
    fprintf(fileID,"BASIC MATERIAL DATA\n\n");
    fprintf(fileID,"MODULUS (Pa): \t\t\t\t%0.6f\n",inputMaterialDataStruct.modulus);
    fprintf(fileID,"COMPRESSIVE YIELD STRESS (Pa): \t\t%0.6f\n",inputMaterialDataStruct.compressiveYieldStress);
    fprintf(fileID,"ULTIMATE COMPRESSIVE STRESS (Pa): \t%0.6f\n",inputMaterialDataStruct.ultimateCompressiveStress);
    fprintf(fileID,"MATERIAL MAXIMUM STRAIN: \t\t%0.6f\n",inputFoamData(end,1));
    fprintf(fileID,"MATERIAL PLASTIC TOUGHNESS ENERGY: \t%0.6f\n",inputMaterialDataStruct.plasticToughnessEnergy);

    fprintf(fileID,"\n\nPLASTIC STRESS STRAIN DATA\n\n");
    fprintf(fileID,"StressData (Pa)\tStrainData\n\n");
    fprintf(fileID,"[\n");
    for(fileWriteIterator = 1:inputDataSize(1))
        for(streamIterator = inputDataSize(2):-1:1)
            switch(streamIterator)
                case 2
                    fprintf(fileID,"[float(%.6f),\t",inputFoamData(fileWriteIterator,streamIterator));
                otherwise
                    fprintf(fileID,"float(%.6f)],\t",inputFoamData(fileWriteIterator,streamIterator));
            end
        end
        fprintf(fileID,"\n");
    end
    fprintf(fileID,"]\n");
    file_close_comparator = fclose(fileID);
    switch(file_close_comparator)
        case(-1)
            fprintf(strcat("File: ",fileName," ", "closing failed\n"));
            return;
        case(0)
            fprintf(strcat("File: ",fileName," ", "closed successfully\n"));
        otherwise
            assert(1==0,"[FILE CLOSE ID SWITCH ERROR]: OUTSIDE CASE STRUCTURE");
    end
end


%Plotting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function success_variable = dartec_unique_plotting(figure_number, title_name, file_name, input_data,x_data_name,x_stream_select, y_data_name, y_stream_select)
    success_variable = -1;
    figure(figure_number)
    hold on
    title_string = strcat(title_name ,"  ", file_name);
    title(sprintf(title_string))
    plot(input_data(:,x_stream_select),input_data(:,y_stream_select))
    xlabel(x_data_name)
    ylabel(y_data_name)
    success_variable = 0;
end
function success_variable = dartec_grouped_plotting(figure_number, title_name, input_data,x_data_name,x_stream_select, y_data_name, y_stream_select,colour_charecter)
    
    figure(figure_number)
    hold on
    title_string = strcat(title_name);
    title(sprintf(title_string))
    plot(input_data(:,x_stream_select),input_data(:,y_stream_select),colour_charecter)
    xlabel(x_data_name)
    ylabel(y_data_name)
    success_variable = 0;
end

%Log Functions
function [] = logDisp(inputText,inputDisplayCondition,DEBUG)
    if(DEBUG)
        if(~isa(inputText,"string"))
            inputText = num2str(inputText);
        end
        fprintf("\n[DEBUG DISPLAY CODE]:%d - %s\n",inputDisplayCondition, inputText);
    end
end
function [] = logWarn(inputText,inputDisplayCondition,DEBUG)
    if(DEBUG)
        if(~isa(inputText,"string"))
            inputText = num2str(inputText);
        end
        fprintf("\n[DEBUG WARNING CODE]:%d - %s\n",inputDisplayCondition, inputText);
    end
end
function [] = logErr(inputText,inputDisplayCondition,DEBUG)
    if(DEBUG)
        if(~isa(inputText,"string"))
            inputText = num2str(inputText);
        end
        fprintf("\n[DEBUG ERROR CODE]:%d - %s\n",inputDisplayCondition, inputText);
    end
end
