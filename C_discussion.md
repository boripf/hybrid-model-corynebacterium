### Mechanistic Model and Parameter Estimation

**Discussion for Batch No. 1**

The mechanistic model for the first batch covers the overall trends of all the parameters. The glucose concentration has the best fit with the rapid glucose consumption during the batch phase and the accumulation during the fed-batch phase. More difficult is the modeling of the biomass and $CO_2$ concentration. With the chosen parameters, the exponential growth during the batch phase is reflected, but the growth stagnation in the following ca. 15-20h is not well covered. While the mechanistic model shows the growth limitation for a few hours, it quickly expects the cells to grow again due to the high concentrations in the end of the fermentation. The modeled concentration of the $CO_2$ follows strongly the trend of the biomass concentration due to its definition within the equations. Concluding, the overall trend is represented by the mechanistic model but it lacks details and is not able to adapt to unexpected and sudden changes due to its constant parameters.

<img src="images\batch_no1_mechanistic.png"
     alt="Batch No.1"
     style="float: center"
     width="600" />

**Discussion of Batch No. 2 with constant parameters througout both phases**

For the mechanistic model with constant parameters throughout both phases, there was a general trade-off to consider when estimating the parameters. During the batch phase, the growth rate and the glucose uptake rate are quite low, while during the fed-batch phase, both parameters increase strongly. Consequently, one of the two phases is well represented. The choice of parameters was based on the good fit of the batch phase, which shows an exponential growth in the beginning and a decrease in growth after all the glucose has been consumed. Therefore, the concentration for both biomass and $CO_2$ are too low in the fed-batch phase. In addition, the lag parameter has the effect of modeling the batch phase well, while it has a large effect on the increase of glucose concentration at the end of the fed-batch phase.

<img src="images\batch_no2_mechanistic.png"
     alt="Batch No.1"
     style="float: center"
     width="600" />

**Discussion of Batch No. 2 with switching parameters for phase 1 & 2**

As a modification of the above mechanistic model, a switch was added to define one parameter set for the batch phase and one for the fed-batch phase. The resulting model reflects the experimental data for biomass, glucose and $CO_2$ concentration very well and could reduce the RMSE from 2.9 to 0.4. Since the $CO_2$ concentration is strongly influenced by the biomass it can still not cover the peak around 9.5h and the decay afterwards.

<img src="images\batch_no2_mechanistic_2stages.png"
     alt="Batch No.1"
     style="float: center"
     width="600" />