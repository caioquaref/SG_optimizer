# GOFAST Automated Correlation Tool
### **Objetive**
I want to develop a correlation tool, that will be able to optimize parameter values to match simulated traces to test data.

---
## Development Notes
- Got GOFAST to run from a python file, YAY!
- Milton Hubbard seems to have good knowledge on optimization programming.
- Wrote a program that can alter any given single parameter. Tables, like electric motor losses, are still a mistery on how to change values in a logical/efficient manner. There is also the fact that some parameters have very few rows and columns, while others have many, this I believe will change the approach to altering them.
- 

---
## To Do List


---
## Backlog
1. Figure out how to treat test data to correlate to the average behavior of the vehicle. This should not be necessary since there isn't always many good test datasets.
2. I need to figure out which parameters I **can** and **can't** change for an optimization. E.g. should I waste time figuring out how to change a parameter that in the end shouldn't be messed with to start?
I should talk to some one with more expirience on the matter, like Nic, Sagar, or maybe Milton.
3. I need to understand how to evaluate the correlation quality between two curves. Refer to paper on references. How to score evaluate and score a simulation result? Will it change depending on maneuver, vehicle architecture, or other?

---
### Quick'n Dirty Backlog
- make GOFAST run
- change files
- iterate runs
- setup optimizer
- select optimization type and define problem
