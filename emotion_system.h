#ifndef EMOTION_SYSTEM_H
#define EMOTION_SYSTEM_H

#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <numeric>

enum class Emotion {
    CALM,      
    EXCITED,   
    SAD,        
    JOYFUL      
};
struct ColorPalette {
    sf::Color primary;
    sf::Color secondary;
};
struct AnimationParameters {
    float speed = 1.0f;
    float intensity = 1.0f;
    float chaos = 0.0f;
};
struct AudioFeatures {
    float volume = 0.0f;
    float timestamp = 0.0f;
    std::vector<float> spectrum; 
    std::vector<float> bandEnergies; 
};

// ---------- 情绪检测器 ----------
class EmotionDetector {
private:
    Emotion currentEmotion = Emotion::CALM;

    float lastVolume = 0.0f;
    std::vector<float> lastSpectrum;

    float smoothVolume = 0.0f;
    float smoothCentroid = 0.0f;
    float smoothFlux = 0.0f;
    float smoothLowHighRatio = 1.0f;
    float smoothLowVar = 0.0f;

    const float alpha = 0.2f;

    float computeSpectralCentroid(const std::vector<float>& spectrum, float sampleRate) {
        float numerator = 0.0f, denominator = 0.0f;
        int N = spectrum.size();
        float binWidth = sampleRate / (2.0f * N);
        for (int i = 0; i < N; ++i) {
            float freq = i * binWidth;
            numerator += freq * spectrum[i];
            denominator += spectrum[i];
        }
        return (denominator > 1e-6) ? numerator / denominator : 0.0f;
    }

    float computeSpectralFlux(const std::vector<float>& current, const std::vector<float>& prev) {
        if (prev.empty() || current.size() != prev.size()) return 0.0f;
        float flux = 0.0f;
        for (size_t i = 0; i < current.size(); ++i) {
            float diff = current[i] - prev[i];
            flux += diff * diff;
        }
        return std::sqrt(flux / current.size());
    }

    float computeLowHighRatio(const std::vector<float>& bandEnergies) {
        if (bandEnergies.size() < 4) return 1.0f;
        int total = bandEnergies.size();
        int lowCount = total / 3;
        int highCount = total / 3;
        float lowSum = 0.0f, highSum = 0.0f;
        for (int i = 0; i < lowCount; ++i) lowSum += bandEnergies[i];
        for (int i = total - highCount; i < total; ++i) highSum += bandEnergies[i];
        float lowAvg = lowSum / lowCount;
        float highAvg = highSum / highCount;
        return (highAvg > 1e-6) ? lowAvg / highAvg : 1.0f;
    }

    float computeLowFreqVariance(const std::vector<float>& bandEnergies) {
        if (bandEnergies.size() < 4) return 0.0f;
        int lowCount = bandEnergies.size() / 3;
        std::vector<float> lowEnergies(lowCount);
        for (int i = 0; i < lowCount; ++i) lowEnergies[i] = bandEnergies[i];
        float mean = std::accumulate(lowEnergies.begin(), lowEnergies.end(), 0.0f) / lowCount;
        float sq_sum = std::inner_product(lowEnergies.begin(), lowEnergies.end(), lowEnergies.begin(), 0.0f);
        return sq_sum / lowCount - mean * mean;
    }

public:
    float getSmoothVolume() const { return smoothVolume; }
    float getSmoothCentroid() const { return smoothCentroid; }
    float getSmoothFlux() const { return smoothFlux; }
    float getSmoothLowHighRatio() const { return smoothLowHighRatio; }
    /*void update(float dt, const AudioFeatures& features, float sampleRate) {
        float volume = features.volume;
        float centroid = computeSpectralCentroid(features.spectrum, sampleRate);
        float flux = computeSpectralFlux(features.spectrum, lastSpectrum);
        float lowHighRatio = computeLowHighRatio(features.bandEnergies);
        float lowVar = computeLowFreqVariance(features.bandEnergies);

        smoothVolume = smoothVolume * (1 - alpha) + volume * alpha;
        smoothCentroid = smoothCentroid * (1 - alpha) + centroid * alpha;
        smoothFlux = smoothFlux * (1 - alpha) + flux * alpha;
        smoothLowHighRatio = smoothLowHighRatio * (1 - alpha) + lowHighRatio * alpha;
        smoothLowVar = smoothLowVar * (1 - alpha) + lowVar * alpha;

        float scoreCalm = 0.0f, scoreExcited = 0.0f, scoreSad = 0.0f, scoreJoyful = 0.0f;

        if (smoothVolume < 0.15f) {
            scoreCalm += 2.0f;
        }
        else if (smoothVolume > 0.3f) {
            if (smoothCentroid > 1000.0f) {
                scoreExcited += 3.0f;
            }
            else {
                scoreExcited += 1.0f;
            }
        }
        else {
            scoreJoyful += 2.0f;
        }

        if (smoothCentroid < 400.0f) {
            scoreSad += 2.0f;
        }
        else if (smoothCentroid > 1400.0f) {
            if (smoothFlux > 0.1f) {
                scoreExcited += 3.0f;
            }
        }
        else {
            scoreJoyful += 2.0f;
        }

        if (smoothFlux < 0.05f) {
            scoreCalm += 2.0f;
        }
        else if (smoothFlux > 0.15f) {
            scoreExcited += 3.0f;
        }
        else {
            scoreJoyful += 2.0f;
        }

        if (smoothLowHighRatio > 1.5f) {
            scoreSad += 2.0f;
        }
        else if (smoothLowHighRatio < 0.8f) {
            scoreExcited += 1.5f;
        }

        if (smoothLowVar < 0.01f) {
            scoreSad += 1.0f;
        }
        else if (smoothLowVar > 0.05f) {
            scoreExcited += 1.5f;
        }

        if (smoothCentroid > 1200.0f && smoothFlux > 0.12f) {
            scoreExcited += 2.0f;
        }
        if (smoothCentroid < 400.0f && smoothFlux < 0.08f) {
            scoreSad += 2.0f;
        }

        if (smoothCentroid > 400.0f && smoothCentroid < 1000.0f &&
            smoothFlux > 0.05f && smoothFlux < 0.15f &&
            smoothLowHighRatio > 1.0f) {
            scoreSad += 3.0f;
        }
        if (smoothLowVar > 0.01f && smoothLowVar < 0.04f) {
            scoreSad += 1.0f;
        }
        if (smoothCentroid > 1000.0f && smoothCentroid < 1400.0f &&
            smoothFlux < 0.05f && smoothLowHighRatio > 1.5f) {
            scoreSad += 3.0f; 
        }

        std::vector<float> scores = { scoreCalm, scoreExcited, scoreSad, scoreJoyful };
        printf("Vol=%.2f, Cent=%.1f, Flux=%.3f, LHRatio=%.2f, LowVar=%.4f\n",
            smoothVolume, smoothCentroid, smoothFlux, smoothLowHighRatio, smoothLowVar);
        printf("Calm=%.1f, Excited=%.1f, Sad=%.1f, Joyful=%.1f -> %s\n",
            scoreCalm, scoreExcited, scoreSad, scoreJoyful, getEmotionName().c_str());
        int idx = std::max_element(scores.begin(), scores.end()) - scores.begin();
        currentEmotion = static_cast<Emotion>(idx);

        lastSpectrum = features.spectrum;
        lastVolume = volume;
    }
    */
    void update(float dt, const AudioFeatures& features, float sampleRate) {
        float volume = features.volume;
        float centroid = computeSpectralCentroid(features.spectrum, sampleRate);
        float flux = computeSpectralFlux(features.spectrum, lastSpectrum);
        float lowHighRatio = computeLowHighRatio(features.bandEnergies);
        float lowVar = computeLowFreqVariance(features.bandEnergies);

        smoothVolume = smoothVolume * (1 - alpha) + volume * alpha;
        smoothCentroid = smoothCentroid * (1 - alpha) + centroid * alpha;
        smoothFlux = smoothFlux * (1 - alpha) + flux * alpha * 10;
        smoothLowHighRatio = smoothLowHighRatio * (1 - alpha) + lowHighRatio * alpha;
        smoothLowVar = smoothLowVar * (1 - alpha) + lowVar * alpha;

        float scoreCalm = 0.0f, scoreExcited = 0.0f, scoreSad = 0.0f, scoreJoyful = 0.0f;

        // ========== 音量权重（Calm 从 2 分降为 1.5） ==========
        if (smoothVolume < 0.15f) {
            scoreCalm += 2.0f;
        }
        else if (smoothVolume > 0.3f) {
            if (smoothCentroid > 1000.0f) {
                scoreExcited += 3.0f;
            }
            else {
                scoreExcited += 1.0f;
            }
        }
        else {
            scoreJoyful += 2.0f;
        }

        // 质心判断
        if (smoothCentroid < 400.0f) {
            scoreSad += 1.0f;
        }
        else if (smoothCentroid > 1400.0f) {
            if (smoothFlux > 0.008f) {
                scoreExcited += 2.0f;
            }
            else if (smoothCentroid > 1600.0f) {
                scoreExcited += 1.0f;
            }
        }
        else {
            scoreJoyful += 2.0f;
        }

        // ========== 谱通量（Calm 从 2 分降为 1.5） ==========
        if (smoothFlux < 0.05f) {
            if (scoreJoyful > 2.0f) scoreCalm += 2.0f;
            else { scoreCalm += 1.0f; }
        }
        else if (smoothFlux > 0.015f) {
            scoreExcited += 3.0f;
        }
        else if (smoothFlux > 0.008f) {
            scoreExcited += 1.0f;
        }
        else {
            scoreJoyful += 2.0f;
        }

        // 低频/高频比
        if (smoothLowHighRatio < 1000.0f) {
            scoreSad += 2.0f;
        }
        else if (smoothLowHighRatio < 0.8f) {
            scoreExcited += 1.5f;
        }

        // 低频稳定性
        if (smoothLowVar < 0.03f) {
            scoreSad += 1.0f;
            if (smoothVolume > 0.05f) scoreSad += 1.0f;
        }
        else if (smoothLowVar > 0.05f) {
            scoreExcited += 1.5f;
        }

        // 特殊组合
        if (smoothCentroid > 1200.0f && smoothFlux > 0.12f) {
            scoreExcited += 2.0f;
        }
        if (smoothCentroid < 400.0f && smoothFlux < 0.08f) {
            scoreSad += 2.0f;
        }

        // ========== 合并后的 Sad 中高音区规则 ==========
        // 条件：质心 400–1400 Hz，通量 < 0.05，低频比 > 1.5 → 加 3 分
        if (smoothCentroid > 400.0f && smoothCentroid < 1400.0f &&
            smoothFlux < 0.05f && smoothLowHighRatio > 1.5f&&scoreExcited<=2.0f) {
            scoreSad +=  1.0f;
        }
        if (scoreExcited >= 2.0f) scoreSad -= 1.0f;

        std::vector<float> scores = { scoreCalm, scoreExcited, scoreSad, scoreJoyful };
        printf("Vol=%.2f, Cent=%.1f, Flux=%.3f, LHRatio=%.2f, LowVar=%.4f\n",
            smoothVolume, smoothCentroid, smoothFlux, smoothLowHighRatio, smoothLowVar);
        printf("Calm=%.1f, Excited=%.1f, Sad=%.1f, Joyful=%.1f -> %s\n",
            scoreCalm, scoreExcited, scoreSad, scoreJoyful, getEmotionName().c_str());
        int idx = std::max_element(scores.begin(), scores.end()) - scores.begin();
        currentEmotion = static_cast<Emotion>(idx);

        lastSpectrum = features.spectrum;
        lastVolume = volume;
    }
    
    Emotion getEmotion() const { return currentEmotion; }

    std::string getEmotionName() const {
        switch (currentEmotion) {
        case Emotion::CALM:    return "Calm";
        case Emotion::EXCITED: return "Excited";
        case Emotion::SAD:     return "Sad";
        case Emotion::JOYFUL:  return "Joyful";
        default: return "Unknown";
        }
    }
};

// ---------- RGB ↔ HSV 转换辅助函数 ----------
namespace {
    void rgbToHsv(sf::Color rgb, float& h, float& s, float& v) {
        float r = rgb.r / 255.0f;
        float g = rgb.g / 255.0f;
        float b = rgb.b / 255.0f;
        float max = std::max({ r, g, b });
        float min = std::min({ r, g, b });
        float delta = max - min;

        v = max;
        if (max != 0.0f)
            s = delta / max;
        else
            s = 0.0f;

        if (delta == 0.0f) {
            h = 0.0f;
        }
        else {
            if (max == r) {
                h = 60.0f * (g - b) / delta;
            }
            else if (max == g) {
                h = 60.0f * (b - r) / delta + 120.0f;
            }
            else {
                h = 60.0f * (r - g) / delta + 240.0f;
            }
            if (h < 0) h += 360.0f;
        }
    }

    sf::Color hsvToRgb(float h, float s, float v) {
        float c = v * s;
        float x = c * (1 - std::abs(std::fmod(h / 60.0f, 2) - 1));
        float m = v - c;

        float r, g, b;
        if (h < 60) { r = c; g = x; b = 0; }
        else if (h < 120) { r = x; g = c; b = 0; }
        else if (h < 180) { r = 0; g = c; b = x; }
        else if (h < 240) { r = 0; g = x; b = c; }
        else if (h < 300) { r = x; g = 0; b = c; }
        else { r = c; g = 0; b = x; }

        return sf::Color(
            static_cast<sf::Uint8>((r + m) * 255),
            static_cast<sf::Uint8>((g + m) * 255),
            static_cast<sf::Uint8>((b + m) * 255)
        );
    }
}

// ---------- 情绪到颜色的映射 ----------
class ColorMapper {
private:
    EmotionDetector& detector;
    float transitionSpeed = 5.0f;
    sf::Color currentPrimary;
    sf::Color currentSecondary;

    sf::Color baseColorFor(Emotion e) const {
        switch (e) {
        case Emotion::CALM:    return sf::Color(150, 210, 255); 
        case Emotion::EXCITED: return sf::Color(255, 150, 100); 
        case Emotion::SAD:     return sf::Color(170, 140, 230); 
        case Emotion::JOYFUL:  return sf::Color(255, 255, 150); 
        default: return sf::Color::White;
        }
    }

    void computeTargetHSV(float& h, float& s, float& v) const {
        Emotion e = detector.getEmotion();
        sf::Color base = baseColorFor(e);

        float baseH, baseS, baseV;
        rgbToHsv(base, baseH, baseS, baseV);

        float vol = detector.getSmoothVolume();
        float cent = detector.getSmoothCentroid();
        float flux = detector.getSmoothFlux();

        float hueOffset = (cent - 1000.0f) * 0.02f;
        hueOffset = std::clamp(hueOffset, -20.0f, 20.0f);
        h = baseH + hueOffset;
        if (h < 0) h += 360.0f;
        else if (h >= 360.0f) h -= 360.0f;

        s = baseS * (0.8f + 0.4f * vol);
        s = std::clamp(s, 0.5f, 1.0f);

        v = baseV * (0.8f + 0.4f * std::min(flux * 3.0f, 1.0f));
        v = std::clamp(v, 0.5f, 1.0f);
    }

    void computeSecondaryHSV(float primaryH, float primaryS, float primaryV,
        float& h2, float& s2, float& v2) const {
        h2 = primaryH + 20.0f;
        if (h2 >= 360.0f) h2 -= 360.0f;
        s2 = primaryS * 0.8f;
        v2 = primaryV * 0.9f;
    }

public:
    ColorMapper(EmotionDetector& det) : detector(det) {
        sf::Color base = baseColorFor(Emotion::CALM);
        currentPrimary = base;
        float h, s, v;
        rgbToHsv(base, h, s, v);
        float h2, s2, v2;
        computeSecondaryHSV(h, s, v, h2, s2, v2);
        currentSecondary = hsvToRgb(h2, s2, v2);
    }

    void update(float dt) {
        float targetH, targetS, targetV;
        computeTargetHSV(targetH, targetS, targetV);

        float targetH2, targetS2, targetV2;
        computeSecondaryHSV(targetH, targetS, targetV, targetH2, targetS2, targetV2);

        sf::Color targetPrimary = hsvToRgb(targetH, targetS, targetV);
        sf::Color targetSecondary = hsvToRgb(targetH2, targetS2, targetV2);

        float t = std::min(1.0f, dt * transitionSpeed);
        currentPrimary.r = currentPrimary.r + (targetPrimary.r - currentPrimary.r) * t;
        currentPrimary.g = currentPrimary.g + (targetPrimary.g - currentPrimary.g) * t;
        currentPrimary.b = currentPrimary.b + (targetPrimary.b - currentPrimary.b) * t;

        currentSecondary.r = currentSecondary.r + (targetSecondary.r - currentSecondary.r) * t;
        currentSecondary.g = currentSecondary.g + (targetSecondary.g - currentSecondary.g) * t;
        currentSecondary.b = currentSecondary.b + (targetSecondary.b - currentSecondary.b) * t;
    }

    ColorPalette getCurrentPalette() const { return { currentPrimary, currentSecondary }; }
    std::string getEmotionName() const { return detector.getEmotionName(); }
};

class AnimationMapper {
private:
    EmotionDetector& detector;
public:
    AnimationMapper(EmotionDetector& det) : detector(det) {}

    AnimationParameters getCurrentParameters() const {
        AnimationParameters params;
        Emotion e = detector.getEmotion();
        switch (e) {
        case Emotion::CALM:
            params.speed = 0.5f; params.intensity = 0.4f; params.chaos = 0.1f; break;
        case Emotion::EXCITED:
            params.speed = 2.0f; params.intensity = 1.8f; params.chaos = 0.7f; break;
        case Emotion::SAD:
            params.speed = 0.3f; params.intensity = 0.6f; params.chaos = 0.2f; break;
        case Emotion::JOYFUL:
            params.speed = 1.5f; params.intensity = 1.3f; params.chaos = 0.5f; break;
        }
        return params;
    }
};
#endif