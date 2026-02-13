#ifndef EMOTION_SYSTEM_H
#define EMOTION_SYSTEM_H

#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <numeric>

// 情绪类型
enum class Emotion {
    CALM,       // 平静
    EXCITED,    // 激昂
    SAD,        // 悲伤
    JOYFUL      // 欢快
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
    std::vector<float> spectrum;          // 完整频谱（用于计算质心、通量）
    std::vector<float> bandEnergies;       // 各频带能量（用于计算高低频比）
};

// ---------- 情绪检测器（基于多特征加权）----------
class EmotionDetector {
private:
    Emotion currentEmotion = Emotion::CALM;

    // 特征缓存（用于计算变化量）
    float lastVolume = 0.0f;
    std::vector<float> lastSpectrum;

    // 平滑后的特征值（用于稳定判断）
    float smoothVolume = 0.0f;
    float smoothCentroid = 0.0f;
    float smoothFlux = 0.0f;
    float smoothLowHighRatio = 1.0f;
    float smoothLowVar = 0.0f;

    const float alpha = 0.2f; // 平滑系数

    // 特征提取函数
    float computeSpectralCentroid(const std::vector<float>& spectrum, float sampleRate) {
        // 质心 = sum(f * magnitude) / sum(magnitude)
        float numerator = 0.0f, denominator = 0.0f;
        int N = spectrum.size();
        float binWidth = sampleRate / (2.0f * N); // 每个bin对应的频率宽度（假设频谱长度为N/2，但这里spectrum已经是N/2个点）
        // 注意：spectrum长度是FFT_SIZE/2，对应频率0~nyquist
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
            flux += diff * diff; // 平方和
        }
        return std::sqrt(flux / current.size()); // 归一化
    }

    float computeLowHighRatio(const std::vector<float>& bandEnergies) {
        if (bandEnergies.size() < 4) return 1.0f;
        int total = bandEnergies.size();
        int lowCount = total / 3;      // 低频区（约0-1/3）
        int highCount = total / 3;      // 高频区（约2/3-末尾）
        float lowSum = 0.0f, highSum = 0.0f;
        for (int i = 0; i < lowCount; ++i) lowSum += bandEnergies[i];
        for (int i = total - highCount; i < total; ++i) highSum += bandEnergies[i];
        float lowAvg = lowSum / lowCount;
        float highAvg = highSum / highCount;
        return (highAvg > 1e-6) ? lowAvg / highAvg : 1.0f; // 低频/高频比
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
    void update(float dt, const AudioFeatures& features, float sampleRate) {
        // 1. 提取当前特征
        float volume = features.volume;
        float centroid = computeSpectralCentroid(features.spectrum, sampleRate);
        float flux = computeSpectralFlux(features.spectrum, lastSpectrum);
        float lowHighRatio = computeLowHighRatio(features.bandEnergies);
        float lowVar = computeLowFreqVariance(features.bandEnergies);

        // 2. 平滑特征（减少毛刺）
        smoothVolume = smoothVolume * (1 - alpha) + volume * alpha;
        smoothCentroid = smoothCentroid * (1 - alpha) + centroid * alpha;
        smoothFlux = smoothFlux * (1 - alpha) + flux * alpha;
        smoothLowHighRatio = smoothLowHighRatio * (1 - alpha) + lowHighRatio * alpha;
        smoothLowVar = smoothLowVar * (1 - alpha) + lowVar * alpha;

        // 3. 根据特征计算各情绪得分（权重需要你根据实际音乐微调）
        float scoreCalm = 0.0f, scoreExcited = 0.0f, scoreSad = 0.0f, scoreJoyful = 0.0f;

        // 音量权重
        if (smoothVolume < 0.15f) scoreCalm += 2.0f;
        else if (smoothVolume > 0.4f) scoreExcited += 2.0f;
        else if (smoothVolume > 0.25f) scoreJoyful += 1.0f;

        // 质心（明亮度）
        if (smoothCentroid < 500.0f) scoreSad += 2.0f;      // 低频主导 → 悲伤
        else if (smoothCentroid > 2000.0f) scoreExcited += 2.0f; // 高频明亮 → 激昂
        else if (smoothCentroid > 1000.0f) scoreJoyful += 1.0f;

        // 谱通量（变化率）
        if (smoothFlux < 0.05f) scoreCalm += 2.0f;
        else if (smoothFlux > 0.2f) scoreExcited += 2.0f;
        else if (smoothFlux > 0.1f) scoreJoyful += 1.0f;

        // 低频/高频比（低频优势 → 悲伤）
        if (smoothLowHighRatio > 1.5f) scoreSad += 2.0f;
        else if (smoothLowHighRatio < 0.7f) scoreExcited += 1.0f;

        // 低频稳定性（低频变化小 → 悲伤，变化大 → 可能激昂）
        if (smoothLowVar < 0.01f) scoreSad += 1.0f;
        else if (smoothLowVar > 0.05f) scoreExcited += 1.0f;

        // 特殊组合：质心高且通量高 → 激昂，质心低且通量低 → 悲伤
        if (smoothCentroid > 1500.0f && smoothFlux > 0.15f) scoreExcited += 2.0f;
        if (smoothCentroid < 800.0f && smoothFlux < 0.08f) scoreSad += 2.0f;

        // 选择最高分情绪
        std::vector<float> scores = { scoreCalm, scoreExcited, scoreSad, scoreJoyful };
        int idx = std::max_element(scores.begin(), scores.end()) - scores.begin();
        currentEmotion = static_cast<Emotion>(idx);

        // 保存当前频谱供下一帧计算通量
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

// ---------- 情绪到颜色的映射（可自定义）----------
class ColorMapper {
private:
    EmotionDetector& detector;
    float transitionSpeed = 5.0f;
    sf::Color currentPrimary;
    sf::Color currentSecondary;

    sf::Color targetPrimaryFor(Emotion e) {
        switch (e) {
        case Emotion::CALM:    return sf::Color(100, 180, 255); // 淡蓝
        case Emotion::EXCITED: return sf::Color(255, 80, 80);   // 亮红
        case Emotion::SAD:     return sf::Color(120, 100, 200); // 紫蓝
        case Emotion::JOYFUL:  return sf::Color(255, 220, 80);  // 金黄
        default: return sf::Color::White;
        }
    }

    sf::Color targetSecondaryFor(Emotion e) {
        switch (e) {
        case Emotion::CALM:    return sf::Color(150, 210, 255); // 更淡蓝
        case Emotion::EXCITED: return sf::Color(255, 150, 100); // 橙红
        case Emotion::SAD:     return sf::Color(170, 140, 230); // 淡紫
        case Emotion::JOYFUL:  return sf::Color(255, 255, 150); // 淡黄
        default: return sf::Color::White;
        }
    }

public:
    ColorMapper(EmotionDetector& det) : detector(det) {
        currentPrimary = targetPrimaryFor(Emotion::CALM);
        currentSecondary = targetSecondaryFor(Emotion::CALM);
    }

    void update(float dt) {
        Emotion e = detector.getEmotion();
        sf::Color targetP = targetPrimaryFor(e);
        sf::Color targetS = targetSecondaryFor(e);

        float t = std::min(1.0f, dt * transitionSpeed);
        currentPrimary.r = currentPrimary.r + (targetP.r - currentPrimary.r) * t;
        currentPrimary.g = currentPrimary.g + (targetP.g - currentPrimary.g) * t;
        currentPrimary.b = currentPrimary.b + (targetP.b - currentPrimary.b) * t;

        currentSecondary.r = currentSecondary.r + (targetS.r - currentSecondary.r) * t;
        currentSecondary.g = currentSecondary.g + (targetS.g - currentSecondary.g) * t;
        currentSecondary.b = currentSecondary.b + (targetS.b - currentSecondary.b) * t;
    }

    ColorPalette getCurrentPalette() const { return { currentPrimary, currentSecondary }; }
    std::string getEmotionName() const { return detector.getEmotionName(); }
};

// ---------- 动画参数映射（可选）----------
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

#endif // EMOTION_SYSTEM_H
