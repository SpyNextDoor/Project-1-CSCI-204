"""! @brief The AudioProcessor package.
"""
import cmath
import math
##
# @file AudioProcessor.py
#
# @brief This package provides the audio processing functions.
#
# @section description_audioprocessor Description
# This package provides the following audio processing functions:
# - AudioProcessor.audio_generate_sine_wave()
#   - It samples a sine wave sound at a given frequency and amplitude.
# - AudioProcessor.audio_generate_square_wave()
#   - It samples a square wave sound at a given frequency and amplitude.
# - AudioProcessor.audio_generate_sawtooth_wave()
#   - It samples a sawtooth wave sound at a given frequency and amplitude.
# - AudioProcessor.audio_generate_complex_wave()
#   - It samples a complex wave sound at a given frequency and amplitude.
# - AudioProcessor.audio_generate_string_wave()
#   - It samples a string wave sound using the Karplus-Strong algorithm at a given frequency and amplitude.
# - AudioProcessor.audio_note_number_to_freq()
#   - It converts the note number to the corresponding wave frequency.
# - AudioProcessor.audio_stereo_gains()
#   - It compute the left and right channel gains given a stereo pan angle.
# - AudioProcessor.audio_multiply_gain()
#   - It multiplies the input audio data by the input gain.
# - AudioProcessor.audio_rise_fall_envelope()
#   - It applies the rise-fall envelope to the input audio data.
# - AudioProcessor.audio_adsr_envelope()
#   - It applies the ADSR envelope to the input audio data.
# - AudioProcessor.audio_stereo_mix_in()
#   - It mixes a mono channel audio data into the stereo audio data.
#
# @section libraries_audioprocessor Libraries/Modules
# - typing (from the standard library)
#   - access to Tuple
# - math (from the standard library)
#   - access to sqrt, pi, sin, cos, exp, and floor
# - DataStructure
#   - access to DataStructure.Array
#
# @section notes_audioprocessor Notes
# - Comments should be Doxygen compatible.
#
# @section todo_audioprocessor TODO
# Your tasks are to implement the following functions: 
# - AudioProcessor.audio_generate_sine_wave()
# - AudioProcessor.audio_generate_square_wave()
# - AudioProcessor.audio_generate_sawtooth_wave()
# - AudioProcessor.audio_generate_complex_wave()
# - AudioProcessor.audio_generate_string_wave()
# - AudioProcessor.audio_note_number_to_freq()
# - AudioProcessor.audio_stereo_gains()
# - AudioProcessor.audio_multiply_gain()
# - AudioProcessor.audio_rise_fall_envelope()
# - AudioProcessor.audio_adsr_envelope()
# - AudioProcessor.audio_stereo_mix_in()
#
# @section author_audioprocessor Author(s)
# - Created by SingChun Lee on 12/24/2023
# - Modified by SingChun Lee on 12/30/2023
#
# Copyright (c) 2023 Bucknell University. All rights reserved.

from typing import Tuple
from math import sqrt, pi, sin, cos, exp, floor
from DataStructure import Array

## The number of attack samples
adsr_attack_samples = 882
## The number of decay samples
adsr_decay_samples = 882
## The number of release samples
adsr_release_samples = 882

# hi its jolien
def audio_generate_sine_wave(audio_data: Array, freq: float, amp: float, samples_per_sec: int) -> None:
    """! This function generates a sine wave audio data.
    
    Given an input frequency **freq** and an amplitude **amp**, this function samples the continuous sine function \f(\sin\f) to the input Array **audio_data**. The sampling rate is determined by the input samples per second **samples_per_sec**. For each sample in **audio_data** (denote it as the **i**-th sample), this function first computes the current time **t** (dividing **i** by **samples_per_sec**). Then, it computes the current angle **theta** (multiplying **t** by **freq** and \f(2\pi\f)). At last, it sets the **i**-th sample **audio_data[i]** to \f(\sin\f)(**theta**) * **amp**.
    
    @param audio_data The audio data.
    
    @param freq The sine wave frequency.
    
    @param amp The sine wave amplitude.
    
    @param samples_per_sec The number of samples per second.
    """
    length_array = len(audio_data)
    for i in range(length_array):
        audio_data[i] = amp * sin(2 * pi * i / samples_per_sec * freq)

def audio_generate_square_wave(audio_data: Array, freq: float, amp: float, samples_per_sec: int) -> None:
    """! This function generates a square wave audio data.
    
    Given an input frequency **freq** and an amplitude **amp**, this function samples the continuous square function \f(\mathrm{square}\f) to the input Array **audio_data**. The sampling rate is determined by the input samples per second **samples_per_sec**. For each sample in **audio_data** (denote it as the **i**-th sample), this function first computes the current time **t** (dividing **i** by **samples_per_sec**). Then, it computes the current angle **theta** (multiplying **t** by **freq** and \f(2\pi\f)). At last, it sets the **i**-th sample **audio_data[i]** to \f(\mathrm{square}\f)(**theta**) * **amp**. Note that, one can compute \f(\mathrm{square}\f)(**theta**) using \f(\sin\f)(**theta**) as follows:
        \f(
        \mathrm{square}(\theta) = \begin{cases}
            1 & \text{if } \sin(\theta) \geq 0 \\
            -1 & \text{otherwise.}
        \end{cases}
        \f)
    
    @param audio_data The audio data.
    
    @param freq The square wave frequency.
    
    @param amp The square wave amplitude.
    
    @param samples_per_sec The number of samples per second.
    """
    array_length = len(audio_data)
    for i in range(array_length):
        multiplier = sin(2 * pi * i/samples_per_sec * freq)
        if multiplier >= 0:
            audio_data[i] = amp * 1
        else:
            audio_data[i] = -amp * 1
    
# Jolien
def audio_generate_sawtooth_wave(audio_data: Array, freq: float, amp: float, samples_per_sec: int) -> None:
    """! This function generates a sawtooth wave audio data.
    
    Given an input frequency **freq** and an amplitude **amp**, this function samples the continuous sawtooth function \f(\mathrm{sawtooth}\f) to the input Array **audio_data**. The sampling rate is determined by the input samples per second **samples_per_sec**. For each sample in **audio_data** (denote it as the **i**-th sample), this function first computes the current time **t** (dividing **i** by **samples_per_sec**). Then, it computes the current cycle **num_cycles** (multiplying **t** by **freq**) and the current sampling position **sample_pos** (the decimal part of **num_cycles**). At last, it sets the **i**-th sample **audio_data[i]** to \f(\mathrm{sawtooth}\f)(**sample_pos**) * **amp**. Note that, one can compute \f(\mathrm{sawtooth}\f)(**sample_pos**) as follows:
        \f(
        \mathrm{sawtooth}(p) = -1 + 2p
        \f)
    
    @param audio_data The audio data.
    
    @param freq The sawtooth wave frequency.
    
    @param amp The sawtooth wave amplitude.
    
    @param samples_per_sec The number of samples per second.
    """
    # p=a(−1+2(tf−⌊tf⌋))
    audio_length = len(audio_data)
    for i in range(audio_length):
        t = i / samples_per_sec
        num_cycles = t * freq
        num_cycles_floor = int(num_cycles)
        decimal = num_cycles - num_cycles_floor
        audio_data[i] = amp * (-1 + 2 * decimal)

def audio_generate_complex_wave(audio_data: Array, freq: float, amp: float, samples_per_sec: int) -> None:
    """! This function generates a complex wave audio data.
    
    Given an input frequency **freq** and an amplitude **amp**, this function samples the continuous complex sine wave function \f(\mathrm{complex}\f) to the input Array **audio_data**. The sampling rate is determined by the input samples per second **samples_per_sec**. For each sample in **audio_data** (denote it as the **i**-th sample), this function first computes the current time **t** (dividing **i** by **samples_per_sec**). Then, it computes the current angle **theta** (multiplying **t** by **freq** and \f(2\pi\f)). At last, it sets the **i**-th sample **audio_data[i]** to \f(\mathrm{complex}\f)(**theta**) * **amp**. In this project, we define the \f(\mathrm{complex}\f)(**theta**) using \f(\mathrm{\sin}\f)(**theta**) as follows:
        \f(
        \mathrm{complex}(\theta) = \frac{\mathrm{sinwave}(\theta)}{\max(|\mathrm{sinwave}|)}
        \f)
    where \f(\mathrm{sinwave}\f)(**theta**) is defined as:
        \f(
        \mathrm{sinwave}(\theta) = \left( (\sin(\theta) + \frac{1}{2}\sin(2\theta) + \frac{1}{4}\sin(3\theta) + \frac{1}{8}\sin(4\theta) + \frac{1}{16}\sin(5\theta) + \frac{1}{32}\sin(6\theta)) \times \exp(-0.0004\theta)\right)^3
        \f)
    and \f(\max(|\mathrm{sinwave}|)\f) is the maximum of the absolute value of all \f(\mathrm{sinwave}(\theta)\f).
    
    @param audio_data The audio data.
    
    @param freq The complex wave frequency.
    
    @param amp The complex wave amplitude.
    
    @param samples_per_sec The number of samples per second.
    """
    length_array = len(audio_data)
    absolute_sin_wave_list = []
    for i in range(length_array):
        theta = 2 * pi * (i/samples_per_sec) * freq
        sin_wave = ((sin(theta) + 0.5 * sin(2*theta) + 0.25 * sin(3*theta) +
                    0.125*sin(4*theta) + 0.0625 * sin(5 * theta) + 0.03125*sin(6 * theta))
                    * exp(-0.0004 * theta))**3
        absolute_vals = abs(sin_wave)
        absolute_sin_wave_list.append(absolute_vals)
        audio_data[i] = amp * sin_wave
    for i in range(length_array):
        audio_data[i] /= max(absolute_sin_wave_list)


# Jolien
def audio_generate_string_wave(audio_data: Array, freq: float, amp: float, samples_per_sec: int) -> None:
    """! This function generates a string wave audio data.
    
    Given an input frequency **freq** and an amplitude **amp**, this function samples the continuous string function \f(\mathrm{string}\f), using the Karplus-Strong algorithm (https://en.wikipedia.org/wiki/Karplus%E2%80%93Strong_string_synthesis), to the input Array **audio_data**. The sampling rate is determined by the input samples per second **samples_per_sec**. In a word, the Karplus-Strong algorithm continuously and repeatedly averages consecutive random wave samples into the output string samples. In this project, this function implements a modified version of the Karplus-Strong algorithm as follows:
    
    - First, instead of using random samples, the function initializes an Array **wave_samples** of size **samples_per_sec** / **freq**, and calls audio_generate_square_wave() to create **wave_samples**.
    - Then, initialize **prev_idx** to the last index of the wave samples and **cur_idx** to 0.
    - For each sample in **audio_data** (denote it as the **i**-th sample), 
      - set the current audio sample (**audio_data[i]**) to the average of the wave samples (**wave_samples**) at **prev_idx** and **cur_idx**.
      - set the wave samples at **cur_idx** to **audio_data[i]**.
      - set **prev_idx** to **cur_idx**.
      - advance **cur_idx** by 1. If it is out of bounds, reset it to 0. Hint: You may handle the out-of-bound situation by using the modulo operator.
    
    @param audio_data The audio data.
    
    @param freq The string wave frequency.
    
    @param amp The string wave amplitude.
    
    @param samples_per_sec The number of samples per second.
    """
    wave_length = int(samples_per_sec / freq)
    audio_length = len(audio_data)
    wave_samples = Array(wave_length)
    audio_generate_square_wave(wave_samples, freq * 100, amp, samples_per_sec)
    prev_idx = wave_length - 1
    cur_idx = 0
    for i in range(audio_length):
        audio_data[i] = (wave_samples[prev_idx] + wave_samples[cur_idx]) / 2
        wave_samples[cur_idx] = audio_data[i]
        prev_idx = cur_idx
        cur_idx = (cur_idx + 1) % wave_length


def audio_note_number_to_freq(note_number: int) -> float:
    """! This function converts the song note number to wave frequency.
    
    Musical note numbers can be converted to their corresponding wave frequencies (Ref: https://www.inspiredacoustics.com/en/MIDI_note_numbers_and_center_frequencies). In short, given a note number **note_number** (denote by \f(n\f)), this function returns the corresponding wave frequency using this formula: \f(440 \times 2^{\frac{n - 69}{12}}\f).
    
    @param note_number The input note number.
    
    @return The converted wave frequency.
    """
    
    freq = 440 * 2 ** ((note_number - 69) / 12)
    return freq

def audio_stereo_gains(angle: float) -> Tuple[float, float]:
    """! This function takes the angle between two stereo sources, computes and returns the two channel gains.
    
    This function computes appropriate gains for the stereo left and right channels, according to a pan angle **angle** specified in radians. It returns a tuple of two gains using the below formula. Given an angle \f(\theta\f), the left \f(L\f) and right \f(L\f) gains are: \f(L = \frac{\sqrt{2}}{2} (\cos(\theta) + \sin(\theta))\f) and \f(R = \frac{\sqrt{2}}{2} (\cos(\theta) - \sin(\theta))\f).
    
    @param angle The pan angle between the two stereo sources.
    
    @return A tuple of the left and right channel gains.
    """
    
    left = sqrt(2)/2*(cos(angle) + sin(angle))
    right = sqrt(2)/2*(cos(angle) - sin(angle))
    stereo_gains = (left, right)
    
    return stereo_gains

def audio_multiply_gain(audio_data: Array, gain: float) -> None:
    """! This function multiplies the audio data by the input gain.
    
    This function multiples the input Array **audio_data** by the input **gain**. For each sample in **audio_data**, it is multiplied by **gain**.
    
    @param audio_data The audio data.
    
    @param gain The input gain.
    """
    
    length_array = len(audio_data)
    for i in range(length_array):
        audio_data[i] *= gain

def audio_rise_fall_envelope(audio_data: Array) -> None:
    """! This functions applies rise/fall envelope to the input audio data.
    
    One technique to make the wave samples sound more natural is to apply an envelope to modify the amplitudes of the wave at different times. This function applies the rise/fall envelope to the input **audio_data**. Let \f(m\f) be the index of the middle sample of **audio_data**. For each sample in **audio_data** (denote it as the **i**-th sample), this function first computes the current gain \f(g(i)\f) by 
        \f(
        g(i) = \begin{cases}
            \frac{i}{m} & \text{if } i <= m \\
            \frac{\text{length of the audio data} - 1 - i}{\text{length of the audio data} - 1 - m} & \text{otherwise.}
        \end{cases}
        \f)
    and then amplifies the current sample (**audio_data[i]**) by \f(g(i)\f).
    
    @param audio_data The audio data.
    """
    num_samples = len(audio_data)
    middle_idx = num_samples // 2

    for i in range(num_samples):
        if i <= middle_idx and middle_idx != 0:
            gain = i / middle_idx
            # audio_data[i] += gain
        else:
            gain = (num_samples - 1 - i) / (num_samples - 1 - middle_idx)
            # audio_data[i] += gain
        audio_data[i] *= gain

def audio_adsr_envelope(audio_data: Array) -> None:
    """! This function applies the ADSR envelope to the input audio data.
    
    One technique to make the wave samples sound more natural is to apply an envelope to modify the amplitudes of the wave at different times. This function applies the attack/decay/substain/release (ADSR) envelope to the input **audio_data**. Ref: https://en.wikipedia.org/wiki/Envelope_(music). The number of attack, decay, and release samples is defined by **adsr_attack_samples**, **adsr_decay_samples**, and **adsr_release_samples**. The first **adsr_attack_samples** samples are in the attack phase, the next **adsr_decay_samples** samples are in the decay phase, and the last **adsr_release_samples** samples are in the release phase. Otherwise, the samples are in the substain phase. In each phase, the gain is computed differently. For each sample in **audio_data** (denote it as the **i**-th sample), this function first computes the current gain \f(g(i)\f) by 
        \f(
        g(i) = \begin{cases}
            \frac{1.2i}{\textbf{adsr_attack_samples}} & \text{if } i < \textbf{adsr_attack_samples} \\
            1 + \frac{0.2(\textbf{adsr_attack_samples} + \textbf{adsr_decay_samples} + i)}{adsr_decay_samples} & \text{if } \textbf{adsr_attack_samples} \leq i < \textbf{adsr_attack_samples} + \textbf{adsr_decay_samples} \\
            1 & \text{if } \textbf{adsr_attack_samples} + \textbf{adsr_decay_samples} \leq i < \text{length of audio data} - \textbf{adsr_release_samples} \\
            \frac{\text{length of audio data} - 1 - i}{\textbf{adsr_release_samples}} & \text{if } \text{length of audio data} - \textbf{adsr_release_samples} \leq i < \text{length of audio data}
        \end{cases}
        \f)
        
    and then amplifies the current sample (**audio_data[i]**) by \f(g(i)\f). Note that if the lengh of **audio_data** is less than **adsr_attack_samples** + **adsr_decay_samples** + **adsr_release_samples**, this function calls audio_rise_fall_envelope() to apply the rise/fall envelope instead.
    
    @param audio_data The audio data. 
    """

    audio_length = len(audio_data)
    if audio_length < adsr_attack_samples + adsr_decay_samples + adsr_release_samples:
        audio_rise_fall_envelope(audio_data)
    else:
        for i in range(audio_length):
            if i < adsr_attack_samples:
                audio_data[i] *= 1.2 * i / adsr_attack_samples
            elif i < adsr_attack_samples + adsr_decay_samples:
                audio_data[i] *= 1.2 - 0.2 * (i - adsr_attack_samples) / adsr_decay_samples
            elif i < audio_length - adsr_release_samples:
                audio_data[i] *= 1
            elif i < audio_length:
                audio_data[i] *= - (i - audio_length) / adsr_release_samples
            else:
                audio_data[i] *= 0
        
# jolien
def audio_stereo_mix_in(stereo_data: Array, channel_data: Array, which_channel: int) -> None:
    """! This function mixes in an audio channel to the stereo audio data.
    
    This function first checks if the number of stereo audio samples (**stereo_data**) is twice the number of input channel data (**channel_data**). If not, it raises an exception. Otherwise, for each sample in **channel_data** (denote it as the **i**-th sample), using the input **which_channel**, this function mixes **channel_data[i]** into **stereo_data** at the position **i** * 2 + **which_channel**, where **which_channel** equal to 0 is the left channel, and 1 is the right channel. When mixing it in, it adds a new value to the existing value.
    
    @param stereo_data The stereo audio data.
    
    @param channel_data The mono channel audio data to mix in.
    
    @param which_channel The channel to mix in.
    """
    
    channel_length = len(channel_data)
    if len(stereo_data) != 2 * channel_length:
        raise Exception
    for i in range(channel_length):
        stereo_data[i * 2 + which_channel] += channel_data[i]
        
