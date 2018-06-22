# Queue-System-Simulation

Discrete event simulation of a data center. This project is inspired by the research work reported in the article Exact analysis of the M/M/k/setup class of Markov chains via recursive renewal reward. The paper studies a dilemma faced by data centre operators on trading off energy consumption and latency. The key question that the paper asks is whether one should turn an idling server off. The advantage of turning off an idling server is that it can save a lot of energy. However, if an off server is needed again, it has to be setup and this adds latency to job processing. This simulation is done to investigate the best setup time and cooling off time set up.

How to use:
1. See test case for reproduce for the input file formats
2. Copy the test case files into the same folder as wrapper.m and simulation.m
3. Adjust the number of test cases to run appropriately in wrapper.m
4. Run wrapper.m
