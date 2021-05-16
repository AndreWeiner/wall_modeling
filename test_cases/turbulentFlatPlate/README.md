# **Steady Flow Over a 2D Flat Plate**

This case is based on the zero pressure gradient flat plate test described
by the Langley Research Center Turbulence Modeling Resource provided at:
```
    https://turbmodels.larc.nasa.gov/flatplate.html
```
The Reynolds number based on the plate length is of Re<sub>L</sub> = 5×10^6

The '*Allrun*' script applies a selection of tubulence models to a range of y+
meshes for which graphs are created that compare the predictions
against the data from Wieghardt.

For further information please visit:
```    
    https://www.openfoam.com/documentation/guides/latest/doc/verification-validation-turbulent-flat-plate-zpg.html
```

This test case is a basic setting for turbulence models due to its simple geometry. By investigating the case, we are able to calculate and plot skin friction '*C<sub>f</sub>*' at the wall of the plate for each y+ which is compared to the analytical solution based on the empirical data by Wieghardt. y+ values are strongly related to the size of meshes. Therefore, we can show how the skin friction changes depending on how the mesh is resolved. It is generally expected that the more the mesh is resolved, then the more accurate the result will be. There are two approaches to deal with the behavior of the fluid at the wall as follows.
1. With Wall Functions
2. Without Wall Functions

Here, the first case will be explained.
## **2D Turbulent Flat Plate Case with Wall Functions**
When we use the case from tutorials in **OpenFOAM**, the related wall functions for '*k*', '*nut*', and '*omega*' are already set in '*bottomWall*' boundary. Each variable has each wall function as follows.
- '*k*'\
'*k*' is called '*turbulence kinetic energy*' that means the kinetic energy per unit mass of the turbulent fluctuations *u'*. As seen below, '*kqRWallFunction*' is used in '*bottomWall*' boundary. '*kqRWallFunction*' is employed for the situation of zero gradient conditions which correspond to Neumann boundary conditions, whereas '*kLowReWallFunction*' gives the boundary conditions of each position of y+.
```
bottomWall
{
    type            kqRWallFunction;
    value           $internalField;
}
```
- '*nut*'\
'*nut*' is called '*turbulence Eddy viscosity*' that is used to solve the closure problem. By using this term, the Reynolds stress can be expressed similar to the typical stress-viscosity relation. Here, '*nutUSpaldingWallFunction*' is used. This function has complicated terms in order to catch the real behavior in the viscous layer as well as the log-law layer. Therefore, it can be used for all area regardless of y+.
```
bottomWall
{
    type            nutUSpaldingWallFunction;
    value           $internalField;
}
```
- '*omega*'\
'*omega*' is called '*specific turbulence dissipation*' that is the rate at which turbulence kinetic energy '*k*' is converted into thermal internal energy per unit volume and time. For '*omega*', '*omegaWallFunction*' is used. This is the only function for '*omega*', and it can be used for all of the regions (viscous and log-law layers) if the blended option is set to '*true*'.
```
bottomWall
{
    type            omegaWallFunction;
    blended         true;
    value           $internalField;
}
```

The result plots are shown in '*PlotCf.ipynb*'. For y+ = 5 and 10, the calculated skin friction values best fit the empirical data. The values also fit well for y+ = 0.05, whereas the values slightly differ from the empirical data for y+ = 1 and 2. For y+ = 30, 50, and 100, the values well fit the empirical data at the back part of the plate (Re<sub>x</sub> > 0.2x10e7), while they definitely cannot catch the behavior of the data at the front part (Re<sub>x</sub> < 0.2x10e7). It implies that the more accurate wall boundary condition can be applied by wall functions than the condition without wall functions. However, there are still some discrepancies between the calculation and the empirical data. A modification of wall functions is needed to reduce these discrepancies.

## **'*Allrun*' Script**

Having mentioned above, the '*Allrun*' script executes the simulation for a various range of y+. Furthermore, a number of code lines are added for several residual limits of **SIMPLE** algorithm and parallel computing. The important parts of the code are explained as follows. Here, '*kOmegaSST*' is used for the simulation. In addition, a list for mesh grading is declared in order to let all the y+ in first cells from boundary layers be the corresponding y+ values in the list.
```
models="
kOmegaSST
"

declare -A grading_vs_yp
#level 5 gradings
grading_vs_yp[0.05]=50000
grading_vs_yp[1]=2200
grading_vs_yp[2]=950
grading_vs_yp[5]=300
grading_vs_yp[10]=130
grading_vs_yp[30]=30
grading_vs_yp[50]=15
grading_vs_yp[100]=5
```

A loop is executed for various models, but there is only one model ('*kOmegaSST*') that will be investigated. Subsequently, the code gives the model information to the dictionary file '*turbulenceProperties*'.

```
for model in $models
do
    echo "Processing model: $model"

    foamDictionary -entry RAS.RASModel -set "$model" \
        constant/turbulenceProperties > /dev/null
```

This loop is for the various y+. Since there are several y+ values to investigate, a separation of data by making each directory corresponding y+ is needed. Afterward, the initial condition is to be copied from the original '*0.model*' (e.g. '*0.kOmegaSST*' or '*0.kEpsilon*') folder when investigating a different y+. Subsequently, the respective grading values are inserted to the '*blockMesh*' file.
```
    for i in "${!grading_vs_yp[@]}"
    do
        yp=$i
        grading=${grading_vs_yp[$yp]}

        # Make directory for file separation
        mkdir yplus_${yp}

        echo "Mesh: y+ $yp"

        # Remove 0.orig folder and copy again from 0.model
        \rm -rf 0.orig
        \cp -rf 0.${model} 0.orig

        sed "s/GRADING/$grading/g" system/blockMeshDict.template > system/blockMeshDict
```

A loop is needed for various residual limits as follows. In the loop, the simulation is actually executed after all the condition is set. When parallel computing codes are used, a number of processor folders are created for each simulations. Hence, these processor folders should be deleted when every simulation starts to execute. Then, the residual limit is added to the '*fvSolution*' dictionary for each step of the loop.
```
        for j in "${!res_value[@]}"
        do
            # Remove all processor folders to decompose the section for every case
            rm -rf processor*
                
            # Add variable for residual limit value
            res=$j
            echo "Residual control: $res"
            # SIMPLE residual control
            sed "s/RESIDUAL/$res/g" system/fvSolution.template > system/fvSolution
```

Now, mesh is to be generated, and then it is decomposed in accordance with the number of processors mentioned in the '*decomposePar*' dictionary file. Subsequently, initial conditions are copied from '*0.orig*' to each processor folder as these folders are deleted after one simulation is finished. Consequently, the '*simpleFoam*' solver with parallel computing is executed here. An option '*-s*' is employed to generate log files with a name '*${model}\_${yp}\_${res}*' (e.g. '*log.simpleFoam.kOmegaSST\_0.05\_1e-5*').
```                
            runApplication -s ${model}_${yp}_${res} blockMesh
            # Decompose mesh for every case
            runApplication -s ${model}_${yp}_${res} decomposePar
            # Distribute 0.orig directory to processor folders
            restore0Dir -processor              
            # Dummy file for post-processing
            touch post.foam

            # Run simpleFoam solver with parallel computing
            runParallel -s ${model}_${yp}_${res} simpleFoam -parallel
```

These lines of code are needed for post-processing. First of all, all the information from a function object '*solverInfo*' is saved and moved to each corresponding y+ folder as '\**.csv'* files. The related function objects will be explained at a later section. Afterward, other post-processing data is saved and moved to separate folders. Table lists should be added only once on the top line of the csv file, and thus echo is also used once before the last loop is started.
```
            # Change solverInfo.dat to solverInfo.csv
            sed 's/#//g' postProcessing/solverInfo/0/solverInfo.dat > postProcessing/solverInfo/0/solverInfo.csv
            cp postProcessing/solverInfo/0/solverInfo.csv yplus_${yp}/solverInfo_${model}_${yp}_${res}.csv
            # Change folder name to preserve each case data
            mv postProcessing postProcessing-${model}_${yp}_${res}                

            UInf=$(foamDictionary -entry internalField 0.orig/U | sed 's/^.*(\s*\([^ ]*\).*/\1/g')
            nuInf=$(foamDictionary -entry nu constant/transportProperties | sed 's/^.*\s\(.*\);/\1/g')

            echo "UInf = $UInf"

            # Echo only once to avoid overwrite and delete '#' in front of ccx
            echo "ccx tau_xx tau_yy tau_zz" > yplus_${yp}/tauw_${model}_${yp}_${res}.csv
            echo "ccx y+" > yplus_${yp}/yplus_${model}_${yp}_${res}.csv
```

This is the final loop in the '*Allrun*' script which iterates over the processors we use. For each processor folder, the related values (e.g. Cx, tau, y+) from the boundary '*bottomWall*' are to be saved. Finally, the values of Cx, tau, and y+ are saved to separate csv files in order to calculate skin friction '*C<sub>f</sub>*'.
```
            # Use an additional loop for gathering all the results from each processor
            for k in $proc_num
            do
                timeDir=$(foamListTimes -processor -latestTime)

                foamDictionary -entry boundaryField.bottomWall.value -value processor$k/$timeDir/Cx | \
                    sed -n '/(/,/)/p' | sed -e 's/[()]//g;/^\s*$/d' > Cx.$$
                foamDictionary -entry boundaryField.bottomWall.value -value processor$k/$timeDir/wallShearStress | \
                    sed -n '/(/,/)/p' | sed -e 's/[()]//g;/^\s*$/d' > tau.$$
                foamDictionary -entry boundaryField.bottomWall.value -value processor$k/$timeDir/yPlus | \
                    sed -n '/(/,/)/p' | sed -e 's/[()]//g;/^\s*$/d' > yp.$$
            
                # Save to csv file instead of dat file
                paste -d ' ' Cx.$$ tau.$$ >> yplus_${yp}/tauw_${model}_${yp}_${res}.csv                
                paste -d ' ' Cx.$$ yp.$$ >> yplus_${yp}/yplus_${model}_${yp}_${res}.csv
                \rm -f Cx.$$ tau.$$ yp.$$                    
            done
```
## **Function Objects in the Turbulent Flat Plate Case**
In **OpenFOAM**, there is a concept called function objects which execute their own missions simultaneously along with the simulation without being compiled before the execution. In order to use the function objects, we need to write '*functions*' in the dictionary file '*controlDict*', and input the objects that we need to use in the simulation. The following function objects are used in this **Turbulent Flat Plate Case**.
```
functions
{
    minMax
    {
        type          fieldMinMax;
        libs          (fieldFunctionObjects);
        writeControl  timeStep; //writeTime;
        fields        (U);
    }
    yPlus
    {
        type            yPlus;
        libs            (fieldFunctionObjects);
        patches         (fixedWall);
        writeControl    writeTime;
    }    
    solverInfo
    {
        type            solverInfo;
        libs            ("libutilityFunctionObjects.so");
        fields          (U p);
        writeResidualFields no;
        writeControl    writeTime;
    }
    #includeFunc "writeCellCentres"
    #includeFunc "wallShearStress"
```
1. '*minMax*' is used to check the convergence of the simulation. If a minimum and a maximum value do not diverge, we are able to decide that this simulation converges and might go to the right way.
2. '*yplus*' is used to save y+ values for designated patches mentioned in the option.
3. '*solverInfo*' is used to save various information for designated fields such as velocity or pressure mentioned in the option. Information of residuals is also saved to the file, and thus we can extract the information from this file.
4. '*writeCellCentres*' saves coordinates of cell centers not to the '*postProcess*' folder, but to each time step folder.
5. '*wallShearStress*' saves maximum and minimum wall shear stress values in the file.