### A_preprocessing

In the pre-processing section, the fermentation data of batch no. 1 and batch no. 2 were prepared for further analysis. During both fermentations, measurements were taken for the biomass and glucose concentration, rate of glucose feed, air flow in the reactor, composition of the offgases, tempertature, pH and dissolved oxygen. Based on data, a lot of conclusions about the fermentation can be conducted. Therefore, the fermentations will be described and discussed in the following.

Both experiments are segmented into a batch phase and a fed-batch phase. However, the fed-batch phases have quite some differences because of different feeding ramps and different start points for the glucose feed. The start of the fed-batch phase starts after 10h for batch No.1 and after 20h for batch No.2.

**Discussion of Batch No.1**
During fermentation, 20 samples were taken to measure biomass and glucose concentration. In the first 10h, which is the batch phase, the microorganisms grow rapidly and consume all the glucose. After that, glucose feeding starts and rapid cell growth is expected. Due to a lack of sampling, it is not entirely clear how the cells behave, but based on the carbon evolution rate (CER), it is expected that limitations will inhibit growth. If the cells are growing normally, the $CO_2$ concentration would be expected to follow the trend of the glucose feed, since rapid growth correlates with higher $CO_2$ production. Instead, the $CO_2$ concentration stagnates between 10h and 20h despite increasing glucose concentration in the medium. Within the 20-25h process duration, trace nutrients were added and the glucose supply was increased to overcome the growth limitation. 
After adding the nutrients, the stirring speed was adjusted to ensure that the nutrients were evenly distributed. The effect of the boost is clearly visible in the biomass concentration, which shows that the cells start growing rapidly again, which continues even after the glucose feed is stopped after 35h. The temperature and the pH are held constant so that no conclusion can be drawn from those parameters.

<img src="images/batch_no1_exp.png"
     alt="Batch No.1"
     style="float: center"
     width="600" />

**Discussion of Batch No.2**
Here again 20 samples were taken to measure biomass and glucose concentration. The batch phase has a duration of 20h, during which the glucose from the medium is quickly completely consumed, which is reflected by the $CO_2$ concentration in the off-gas. Based on these measurements, *Corynebacterium glutamicum* grows exponentially in the beginning, but stops rapidly after about 9.5h, where the $CO_2$ concentration drops to nearly zero with time. After 20h of process time, glucose feeding is started, which increases cell growth again. The growth rate is so high that the glucose concentration stays at zero for the rest of the fermentation, leading to a high glucose uptake rate. There is no indication in the measurements that growth limitations as in the first fermentation occur in this experiment.

<img src="images/batch_no2_exp.png"
     alt="Batch No.2"
     style="float: center"
     width="600" />

**Conclusion**
- during batch phase, the glucose was quickly consumed in both fermentation experiments
- glucose feed enhances the growth of *Corynebacterium glutamicum*
- growth limitations as in batch #1 were overcome in batch #2
- many parameters are controlled such as pH, T, DO2 leading to limited number of parameters for modeling
- only parameters we can model are biomass, substrate, co2
- for the volume we don't have data to compare it with

Given that our model is not a realistic model, the mechanistic and data-driven models are still defined for the learning experience of developing a hybrid model.
