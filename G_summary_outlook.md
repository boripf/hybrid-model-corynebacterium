### G_Hybrid semi-parametric model

By combining data-driven flexibility with mechanistic model interpretability, the hybrid semi-parametric model is capable of accurate one-step-ahead prediction. This characteristic is particularly true for average batches. When predicting deviations from average data, however, its limitations become apparent. One reason for the limited accuracy is that the training data generated by the Monte Carlo method contains little variability. This is due to the fact that there are only two input parameters. Notably, the glucose feed and the total feed to the reactor are the same for all batches, as they cannot be modeled. To improve model robustness, data generation should include additional parameters that introduce more changes that better reflect real-world scenarios. This is especially important for data-driven models, which thrive on diverse inputs to improve accuracy by adapting to changes.

<img src="images/batch_no1_exp.png"
     alt="Batch No.1"
     style="float: center"
     width="600" />

To further enhance the effectiveness of the model, sensitivity analysis is identified as a key area for refinement. Performing such an analysis on Batch #2 promises deeper insights into influencing factors, shedding light on model dynamics and suggesting avenues for improvement. Understanding the nuances of sensitivity analysis can facilitate more informed parameter adjustments, potentially leading to optimized reaction kinetics. In particular, the insignificance of the lag parameter in the fed-batch simulations for Batch #2 highlights its potential irrelevance, which may warrant a re-evaluation of its inclusion.

Going forward, the evolution of the model could include two important avenues: first, diving deeper into the correlation and identificationability of parameters to provide a more nuanced understanding of their relationships and their role in shaping model results. Parameter identificationability refers to the extent to which the values of model parameters can be accurately estimated or determined from experimental data. It therefore assesses the feasibility of finding parameter values that best match the model predictions with observed experimental results, which is influenced by data quality, model complexity, and parameter uniqueness. Second, iterating through all the development steps with the knowledge gained provides an opportunity to refine the accuracy and predictive power of the model.