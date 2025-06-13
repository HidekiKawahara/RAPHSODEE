# RAPHSODEE

A new series of acoustic attributes analysis. This is an extension of safeguarding sounds and CAPRICEP. This also has a fast but ecologically inefficient real-time reverberation application.

## New important implementation (2025 13 June)
I uploaded a MATLAB m-function to make any signal into a safeguarded signal. Thanks to the giant-FFT SRC, you can make the signal sampling frequency as you like. Please use the interactive and real-time acoustic condition checker tool for testing the generated safeguarded signals. The name of the function is:

signalSafeguardwithGiantFFTSRC

Have fun!

Please check samplingRateConvByDFTwin.m function in src. It is an independent implementation of the following algorithm.
% Välimäki, V., & Bilbao, S. (2023). Giant FFTs for Sample-Rate Conversion.
% Journal of the Audio Engineering Society, 71(3), 88–99.
% DOI: https://doi.org/10.17743/jaes.2022.0061

## Note: retrospective use of safeguarding (2025 June 13)

You don't need to safeguard a signal. You can enjoy benefits of safeguarding. Plase refer the following script in src:

testScrptForEA2025July.m

## (Pre) Release note
### Updated an interactive and real-time tool for testing acoustic conditions (Annual spring meeting of ASJ 2025)

Name of application: acousticConditionChecker

### Interactive and real-time tools for testing acoustic environment effects on speech production: prerelease version (ASJ Hearing Research technical meeting May 2024)

Name of application: reverbeLmbrdSpeechTesterRev

Quick start guide for this tool:
<[https://www.youtube.com/watch?v=K-z5juxZg0o](https://www.youtube.com/watch?v=K-z5juxZg0o)>

### Interactive and real-time tool for acoustic environment analysis: prerelease version, 02/Nov./2023, APSIPA2023 version

Name of application: acousticConditionChecker

Quick start guide for this tool:
<[https://youtu.be/3Ial0qaEBfs?si=H3UbsTpC37mwZRxc](https://youtu.be/3Ial0qaEBfs?si=H3UbsTpC37mwZRxc)>

How to calibrate using a sound pressure meter
<https://youtu.be/s1oBoDsT5zA?si=CISF4i3SQLgEJQ94>

## References

Kawahara, H., Yatabe, K., Matsui, T., Hodoshima, N., Mizumachi, M., and Sakakibara, K. (2024) Real-time processing and interface design of acoustic environment control tool for voice production research, Proc. Auditory Res. Meeting, The Acoustical Society of Japan, Vol.54, No.3, H-2024-48, PP.236-252. (In Japanese).

Kawahara, H., Yatabe, K., Sakakibara, K.-I., Mizumachi, M., & Kitamura, T. (2023). Simultaneous Measurement of Multiple Acoustic Attributes Using Structured Periodic Test Signals, Including Music and Other Sound Materials. In arXiv [cs.SD]. arXiv. <http://arxiv.org/abs/2309.02767>

This article was presented at APSIPA2023 on 1st Nov. 2023. It will be accessible on the APSIPA website.

<http://www.apsipa.org/>
