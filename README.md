# Source_Tract_Coupling
This program is written for c-language. It includes

SourceTractModel.c : Simulate vocal fold oscillations coupled to vocal tract acoustics.

SubglottalAreaF.txt : Area function of subglottis

Compile: gcc SourceTractModel.c -o SourceTractModel -lm

Executable file: SourceTractModel

Output: Time vs glottal airflow is written int0 the output file "Flow.txt" 

Running this code will reproduce partial results in the publication :
I. T. Tokuda, M. Zemke, M. Kob, and H. Herzel, "Biomechanical modeling of register transitions and the role of vocal tract resonators," Journal of the Acoustical Society of America, Vol. 127, Issue 3, pp. 1528-1536 (2010); C. T. Herbst, C. P. H. Elemans, I. T. Tokuda, V. Chatziioannous, J. G. Svec, "Dynamic system coupling in voice production," Frontiers in Network Physiology (2022)

The variable "T" in the main function represents the tension parameter that controls the fundamental frequency. Setting it to "T=3.3" / "T=4.0" will reproduce Fig. 5B / Fig. 5C in Herbst et al. (2022).
