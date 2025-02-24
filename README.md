# Estimation of the Average Treatment Effect (ATE): Practical Recommendations

Authors:

- [Charlotte Voinot](https://chvoinot.github.io/) - [INRIA PreMeDICaL team, Idesp, Université de Montpellier](https://team.inria.fr/premedical/) and [Sanofi R&D](https://www.sanofi.com/en)
- [Clément Berenfeld](https://cberenfeld.github.io/) - [Universität Potsdam, Potsdam, Germany](https://www.uni-potsdam.de/en/university-of-potsdam)
- Imke Mayer - [Charité – Universität Berlin, Berlin, Germany](https://www.charite.de/)
- Bernard Sebastien - [Sanofi R&D](https://www.sanofi.com/en)
- [Julie Josse](http://juliejosse.com/) - [INRIA PreMeDICaL team, Idesp, Université de Montpellier](https://team.inria.fr/premedical/)

Causal survival analysis combines survival analysis and causal inference to evaluate the effect of a treatment or intervention on a time-to-event outcome, such as survival time. It offers an alternative to relying solely on Cox models for assessing these effects. In this paper, we present a comprehensive review of estimators for the average treatment effect measured with the restricted mean survival time, including regression-based methods, weighting approaches, and hybrid techniques. We investigate their theoretical properties and compare their performance through extensive numerical experiments. Our analysis focuses on the finite-sample behavior of these estimators, the influence of nuisance parameter selection, and their robustness and stability under model misspecification. By bridging theoretical insights with practical evaluation, we aim to equip practitioners with both state-of-the-art implementations of these methods and practical guidelines for selecting appropriate estimators for treatment effect estimation. Among the approaches considered, G-formula two-learners, AIPCW-AIPTW, Buckley-James estimators, and causal survival forests emerge as particularly promising.


[![build and publish](https://github.com/computorg/template-computo-R/actions/workflows/build.yml/badge.svg)](https://github.com/computorg/template-computo-R/actions/workflows/build.yml)
