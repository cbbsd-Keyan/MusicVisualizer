#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <complex>
#include "SmoothValue.h"
#include "emotion_system.h"

EmotionDetector emotionDetector;
ColorMapper colorMapper(emotionDetector);
AnimationMapper animationMapper(emotionDetector);

// ==================== 窗口常量 ====================
const int WINDOW_WIDTH = 1200;
const int WINDOW_HEIGHT = 800;

// ==================== 公共工具函数 ====================
float calculateVolume(const sf::Int16* samples, size_t count) {
    if (count == 0) return 0.0f;
    float sumSquares = 0.0f;
    for (size_t i = 0; i < count; i++) {
        float sample = samples[i] / 32768.0f;
        sumSquares += sample * sample;
    }
    return std::sqrt(sumSquares / count);
}

// ---------- 以下是从 circle.cpp 原封不动搬来的 HSB 和 FFT ----------
sf::Color hsbToRgb(float hue, float saturation, float brightness) {
    hue = fmod(hue, 360.0f);
    if (hue < 0) hue += 360.0f;
    saturation = std::clamp(saturation, 0.0f, 100.0f) / 100.0f;
    brightness = std::clamp(brightness, 0.0f, 100.0f) / 100.0f;

    float c = brightness * saturation;
    float x = c * (1.0f - fabs(fmod(hue / 60.0f, 2.0f) - 1.0f));
    float m = brightness - c;

    float r = 0.0f, g = 0.0f, b = 0.0f;

    if (hue >= 0 && hue < 60) {
        r = c; g = x; b = 0;
    }
    else if (hue >= 60 && hue < 120) {
        r = x; g = c; b = 0;
    }
    else if (hue >= 120 && hue < 180) {
        r = 0; g = c; b = x;
    }
    else if (hue >= 180 && hue < 240) {
        r = 0; g = x; b = c;
    }
    else if (hue >= 240 && hue < 300) {
        r = x; g = 0; b = c;
    }
    else {
        r = c; g = 0; b = x;
    }

    return sf::Color(
        static_cast<sf::Uint8>((r + m) * 255),
        static_cast<sf::Uint8>((g + m) * 255),
        static_cast<sf::Uint8>((b + m) * 255)
    );
}

void computeFFT(const std::vector<float>& input, std::vector<float>& spectrum) {
    int N = (int)input.size();
    if (N == 0) return;
    std::vector<std::complex<float>> data(N);
    for (int i = 0; i < N; i++) data[i] = std::complex<float>(input[i], 0.0f);
    int j = 0;
    for (int i = 0; i < N; i++) {
        if (j > i) std::swap(data[i], data[j]);
        int m = N >> 1;
        while (m >= 1 && j >= m) { j -= m; m >>= 1; }
        j += m;
    }
    for (int s = 1; s <= (int)std::log2(N); s++) {
        int m = 1 << s;
        float theta = -2.0f * 3.14159265358979323846f / m;
        std::complex<float> wm(cos(theta), sin(theta));
        for (int k = 0; k < N; k += m) {
            std::complex<float> w(1.0f, 0.0f);
            for (int j = 0; j < m / 2; j++) {
                std::complex<float> t = w * data[k + j + m / 2];
                std::complex<float> u = data[k + j];
                data[k + j] = u + t;
                data[k + j + m / 2] = u - t;
                w *= wm;
            }
        }
    }
    spectrum.resize(N / 2);
    for (int i = 0; i < N / 2; i++) spectrum[i] = std::abs(data[i]) / (N / 2);
}

std::vector<float> mapSpectrumToBands(const std::vector<float>& spectrum, int bands, float sampleRate) {
    std::vector<float> bandEnergies(bands, 0.0f);
    int specSize = (int)spectrum.size();
    if (specSize == 0 || bands == 0) return bandEnergies;

    // ---------- 使用传入的真实采样率 ----------
    float nyquist = sampleRate / 2.0f;
    float minFreq = 20.0f;
    float maxFreq = nyquist;          // 直接使用 nyquist，不硬编码 20000

    float logMin = std::log10(minFreq);
    float logMax = std::log10(maxFreq);

    // ---------- 频带分配（保证不重叠，至少1个bin）----------
    std::vector<int> binStart(bands), binEnd(bands);
    int lastEnd = 0;

    for (int b = 0; b < bands; b++) {
        // 对数均匀的频率范围
        float logStart = logMin + (logMax - logMin) * b / bands;
        float logEnd = logMin + (logMax - logMin) * (b + 1) / bands;
        float freqStart = std::pow(10.0f, logStart);
        float freqEnd = std::pow(10.0f, logEnd);

        // 转换为 bin 索引（浮点）
        float binStartF = freqStart / nyquist * specSize;
        float binEndF = freqEnd / nyquist * specSize;

        // 整数分配，保证不重叠且至少1个bin
        binStart[b] = std::max((int)std::round(binStartF), lastEnd);
        if (b == bands - 1) {
            binEnd[b] = specSize;   // 最后一个频带到结尾
        }
        else {
            binEnd[b] = std::max(binStart[b] + 1, std::min((int)std::round(binEndF), specSize));
        }
        binEnd[b] = std::min(binEnd[b], specSize);
        if (binEnd[b] <= binStart[b]) binEnd[b] = binStart[b] + 1;

        lastEnd = binEnd[b];
    }

    // ---------- 计算每个频带的平均能量 ----------
    std::vector<float> rawEnergies(bands, 0.0f);
    for (int b = 0; b < bands; b++) {
        float sum = 0.0f;
        int count = binEnd[b] - binStart[b];
        for (int i = binStart[b]; i < binEnd[b]; i++) {
            sum += spectrum[i];
        }
        float avg = sum / count;
        rawEnergies[b] = std::log10(1.0f + avg * 99.0f);  // 对数压缩
    }

    // ---------- 自动增益控制（AGC）----------
    static std::vector<float> longTermAvg(bands, 0.0f);
    static std::vector<float> gain(bands, 1.0f);
    float alpha = 0.02f;
    float targetLevel = 0.4f;        // 可调，整体高度

    for (int b = 0; b < bands; b++) {
        longTermAvg[b] = longTermAvg[b] * (1 - alpha) + rawEnergies[b] * alpha;
        if (longTermAvg[b] > 0.01f) {
            gain[b] = targetLevel / longTermAvg[b];
        }
        // 🔥 关键修复：大幅提高增益上限，让高频能跳起来
        gain[b] = std::clamp(gain[b], 0.2f, 8.0f);   // 上限从4.0 → 8.0

        // 🔥 高频额外增益：最后1/10频带再提升
        float extraGain = 1.0f;
        if (b > bands * 0.9f) extraGain = 4.0f;      // 最后10%频带增益4倍
        else if (b > bands * 0.7f) extraGain = 2.5f; // 70%~90%频带增益2.5倍

        bandEnergies[b] = rawEnergies[b] * gain[b] * extraGain;
    }

    return bandEnergies;
}

// ==================== 效果 1：粒子系统（来自 particles.cpp，原封不动）====================
class MusicParticle {
public:
    sf::Vector2f position;
    sf::Vector2f velocity;
    sf::Color color;
    float radius;
    float life;
    MusicParticle(float x, float y) : position(x, y), velocity(0, 0), color(sf::Color::White), radius(10.0f), life(1.0f) {}
    void update(float dt, float volume) {
        life -= dt * 0.5f;
        position += velocity * dt * 60.0f;
        velocity.y += 0.03f;
        radius = 8.0f + volume * 50.0f;
    }
    void draw(sf::RenderTarget& target) {
        if (life <= 0) return;
        sf::CircleShape circle(radius);
        circle.setPosition(position);
        circle.setFillColor(sf::Color(color.r, color.g, color.b, static_cast<sf::Uint8>(life * 255)));
        target.draw(circle);
    }
};

namespace ParticleEffect {
    std::vector<MusicParticle> particles;
    sf::Clock particleClock;
    SmoothValue<float> smoothedVolume(0.0f, 8.0f);
    void init() { particles.clear(); particleClock.restart(); }
    void update(float dt, float volume, bool isPlaying, const ColorPalette& palette, const AnimationParameters& animParams,
        bool useEmotionalColors) {   // 新增参数
        smoothedVolume.setTarget(volume);
        smoothedVolume.update(dt);
        float vol = smoothedVolume.getCurrent();
        if (isPlaying) {
            if (particleClock.getElapsedTime().asSeconds() > 0.05f) {
                int particleCount = 2 + static_cast<int>(vol * 5.0f);
                for (int i = 0; i < particleCount; i++) {
                    float angle = static_cast<float>(rand() % 360) * 3.14159f / 180.0f;
                    float circleRadius = 20.0f + vol * 50.0f;
                    MusicParticle particle(
                        WINDOW_WIDTH / 2.0f + cos(angle) * circleRadius,
                        WINDOW_HEIGHT / 2.0f + sin(angle) * circleRadius
                    );
                    // ---------- 颜色选择 ----------
                    if (useEmotionalColors) {
                        // 情感模式：使用调色板主色
                        particle.color = palette.primary;
                    }
                    else {
                        // 彩色模式：原来的 HSB 颜色（根据角度）
                        float hue = fmod(angle * 180.0f / 3.14159f, 360.0f);
                        particle.color = sf::Color(
                            static_cast<sf::Uint8>(127 + 127 * sin(hue * 3.14159f / 180.0f)),
                            static_cast<sf::Uint8>(127 + 127 * sin((hue + 120.0f) * 3.14159f / 180.0f)),
                            static_cast<sf::Uint8>(127 + 127 * sin((hue + 240.0f) * 3.14159f / 180.0f))
                        );
                    }
                    // 速度（情感模式可乘 animParams.speed，彩色模式保持原样）
                    particle.velocity = sf::Vector2f(
                        cos(angle) * (2.5f + vol * 10.0f) * (useEmotionalColors ? animParams.speed : 1.0f),
                        sin(angle) * (2.5f + vol * 10.0f) * (useEmotionalColors ? animParams.speed : 1.0f)
                    );
                    particles.push_back(particle);
                }
                particleClock.restart();
            }
        }
        for (auto it = particles.begin(); it != particles.end();) {
            it->update(dt, vol);
            if (it->life <= 0) it = particles.erase(it);
            else ++it;
        }
    }
    void draw(sf::RenderWindow& window) { for (auto& p : particles) p.draw(window); }
    void clear() { particles.clear(); }
}

// ==================== 效果 2：波形图（来自 waveform_visualizer_english.cpp，原封不动）====================
class WaveformVisualizer {
private:
    sf::VertexArray waveform;
    std::vector<float> audioBuffer;
    float scaleFactor;
public:
    WaveformVisualizer() : waveform(sf::LineStrip, 500), scaleFactor(100.0f) {
        for (int i = 0; i < 500; i++) {
            float x = static_cast<float>(i) / (500 - 1) * WINDOW_WIDTH;
            waveform[i].position = sf::Vector2f(x, WINDOW_HEIGHT / 2);
            waveform[i].color = sf::Color::White;
        }
    }
    void setColorsFromPalette(const ColorPalette& palette, float volume) {
        for (int i = 0; i < 500; i++) {
            float ratio = (float)i / 500;
            sf::Color base;
            if (ratio < 0.5f) base = palette.primary;
            else base = palette.secondary;
            waveform[i].color = sf::Color(base.r, base.g, base.b, static_cast<sf::Uint8>(150 + volume * 105));
        }
    }
    void update(const std::vector<float>& samples, float volume, float time) {
        if (samples.empty()) return;
        float amplitude = scaleFactor * (5.0f + volume * 15.0f);
        float timeOffset = sin(time * 2.0f) * 10.0f;
        for (int i = 0; i < 500; i++) {
            int sampleIndex = (i * (int)samples.size()) / 500;
            sampleIndex = std::min(sampleIndex, (int)samples.size() - 1);
            float sampleValue = samples[sampleIndex];
            static std::vector<float> prevValues(500, 0.0f);
            float smoothedValue = prevValues[i] * 0.7f + sampleValue * 0.3f;
            prevValues[i] = smoothedValue;
            float x = static_cast<float>(i) / (500 - 1) * WINDOW_WIDTH;
            float y = WINDOW_HEIGHT / 2 + smoothedValue * amplitude + sin(i * 0.1f + time * 3.0f) * (5.0f + volume * 15.0f);
            waveform[i].position = sf::Vector2f(x, y);
            float hue = fmod(i * 0.5f + time * 50.0f, 360.0f);
            float ratio = hue / 60.0f;
            int sector = static_cast<int>(ratio) % 6;
            float fraction = ratio - sector;
            sf::Uint8 r, g, b;
            switch (sector) {
            case 0: r = 255; g = static_cast<sf::Uint8>(fraction * 255); b = 0; break;
            case 1: r = static_cast<sf::Uint8>((1 - fraction) * 255); g = 255; b = 0; break;
            case 2: r = 0; g = 255; b = static_cast<sf::Uint8>(fraction * 255); break;
            case 3: r = 0; g = static_cast<sf::Uint8>((1 - fraction) * 255); b = 255; break;
            case 4: r = static_cast<sf::Uint8>(fraction * 255); g = 0; b = 255; break;
            default: r = 255; g = 0; b = static_cast<sf::Uint8>((1 - fraction) * 255); break;
            }
            sf::Uint8 alpha = static_cast<sf::Uint8>(150 + volume * 105);
            waveform[i].color = sf::Color(r, g, b, alpha);
        }
    }
    void draw(sf::RenderTarget& target) {
        for (int offset = -10; offset <= 10; offset++) {
            sf::VertexArray thickLine = waveform;
            for (unsigned int i = 0; i < thickLine.getVertexCount(); i++) thickLine[i].position.y += static_cast<float>(offset) * 0.3f;
            target.draw(thickLine);
        }
        sf::VertexArray background(sf::Quads, 4);
        background[0].position = sf::Vector2f(0, WINDOW_HEIGHT / 2 - 150);
        background[1].position = sf::Vector2f(WINDOW_WIDTH, WINDOW_HEIGHT / 2 - 150);
        background[2].position = sf::Vector2f(WINDOW_WIDTH, WINDOW_HEIGHT / 2 + 150);
        background[3].position = sf::Vector2f(0, WINDOW_HEIGHT / 2 + 150);
        background[0].color = sf::Color(10, 10, 40, 50);
        background[1].color = sf::Color(10, 10, 40, 50);
        background[2].color = sf::Color(10, 10, 40, 0);
        background[3].color = sf::Color(10, 10, 40, 0);
        target.draw(background);
    }
    void setScaleFactor(float scale) { scaleFactor = scale; }
};

namespace WaveformEffect {
    WaveformVisualizer waveform;
    SmoothValue<float> smoothedVolume(0.0f, 10.0f);
    float currentScale = 100.0f;
    bool colorMode = true;
    std::vector<sf::CircleShape> backgroundParticles;
    void init() {
        for (int i = 0; i < 100; i++) {
            sf::CircleShape particle(1.0f + (rand() % 100) / 100.0f * 3.0f);
            particle.setPosition(rand() % WINDOW_WIDTH, rand() % WINDOW_HEIGHT);
            particle.setFillColor(sf::Color(50, 100, 200, 30));
            backgroundParticles.push_back(particle);
        }
        waveform.setScaleFactor(currentScale);
    }
    void update(float dt, float volume, float time, const std::vector<float>& samples,
        bool isPlaying, const ColorPalette& palette, const AnimationParameters& animParams,
        bool useEmotionalColors) {
        smoothedVolume.setTarget(volume);
        smoothedVolume.update(dt);
        static std::vector<float> decaySamples;
        if (isPlaying) {
            if (!samples.empty()) {
                waveform.update(samples, smoothedVolume.getCurrent(), time);
                decaySamples = samples;
            }
        }
        else {
            if (decaySamples.empty()) {
                decaySamples.assign(500, 0.0f);
            }
            else {
                for (float& val : decaySamples) val *= 0.95f;
            }
            waveform.update(decaySamples, smoothedVolume.getCurrent() * 0.3f, time);
        }
        // ---------- 根据模式设置波形颜色 ----------
        if (useEmotionalColors) {
            waveform.setColorsFromPalette(palette, smoothedVolume.getCurrent());
        }
        else {
        }
        for (auto& particle : backgroundParticles) {
            particle.move(0.1f, 0.05f);
            if (particle.getPosition().y > WINDOW_HEIGHT)
                particle.setPosition(rand() % WINDOW_WIDTH, -10.0f);
        }
    }
    void draw(sf::RenderWindow& window) {
        for (auto& p : backgroundParticles) window.draw(p);
        waveform.draw(window);
        sf::RectangleShape centerLine(sf::Vector2f(WINDOW_WIDTH, 1));
        centerLine.setPosition(0, WINDOW_HEIGHT / 2);
        centerLine.setFillColor(sf::Color(255, 255, 255, 50));
        window.draw(centerLine);
    }
    void setScale(float scale) { currentScale = scale; waveform.setScaleFactor(scale); }
}

// ==================== 效果 3：中央圆环 ====================
namespace CircleEffect {
    sf::CircleShape centerCircle(0.0f);
    SmoothValue<float> circleSize(0.0f, 10.0f);
    SmoothValue<float> smoothedVolume(0.0f, 5.0f);
    sf::RenderTexture trailTexture;
    void init() {
        centerCircle.setOrigin(0, 0);
        centerCircle.setPosition(WINDOW_WIDTH / 2.0f, WINDOW_HEIGHT / 2.0f);
        centerCircle.setFillColor(sf::Color::Transparent);
        centerCircle.setOutlineThickness(3.0f);
        centerCircle.setOutlineColor(sf::Color::White);
        if (!trailTexture.create(WINDOW_WIDTH, WINDOW_HEIGHT)) {
            std::cerr << "Circle: Cannot create trail texture!" << std::endl;
        }
        trailTexture.clear(sf::Color::Black);
    }
    void update(float dt, float volume, const ColorPalette& palette, bool useEmotionalColors) {
        smoothedVolume.setTarget(volume);
        smoothedVolume.update(dt);
        float targetCircleSize = smoothedVolume.getCurrent() * 2200.0f;
        targetCircleSize = std::min(targetCircleSize, 400.0f);
        circleSize.setTarget(targetCircleSize);
        circleSize.update(dt);
        float currentCircleSize = circleSize.getCurrent();
        centerCircle.setRadius(currentCircleSize);
        centerCircle.setOrigin(currentCircleSize, currentCircleSize);

        // ---------- 设置颜色 ----------
        if (useEmotionalColors) {
            centerCircle.setOutlineColor(palette.primary); 
            centerCircle.setFillColor(sf::Color::Transparent); 
        }
        else {
            centerCircle.setOutlineColor(sf::Color::White);
            centerCircle.setFillColor(sf::Color::Transparent);
        }
        sf::RectangleShape trailRect(sf::Vector2f(WINDOW_WIDTH, WINDOW_HEIGHT));
        trailRect.setFillColor(sf::Color(0, 0, 0, 30));
        trailTexture.draw(trailRect, sf::BlendAlpha);
        if (currentCircleSize > 5.0f) {
            trailTexture.draw(centerCircle);
        }
        trailTexture.display();
    }
    void draw(sf::RenderWindow& window) {
        sf::Sprite trailSprite(trailTexture.getTexture());
        window.draw(trailSprite);
    }
}

// ==================== 效果 4：频谱条 ====================
class FrequencyBars {
private:
    std::vector<sf::RectangleShape> bars;
    std::vector<SmoothValue<float>> barHeights;
    int barCount;
    float barWidth;
    float maxHeight;
public:
    FrequencyBars(int count, float maxH = 400.0f) : barCount(count), maxHeight(maxH) {
        float totalSpacing = (count - 1) * 1.0f;
        float availableWidth = WINDOW_WIDTH - totalSpacing;
        barWidth = availableWidth / count;
        for (int i = 0; i < count; i++) {
            sf::RectangleShape bar(sf::Vector2f(barWidth, 10));
            float xPos = i * (barWidth + 1.0f);
            bar.setPosition(xPos, WINDOW_HEIGHT);
            bar.setFillColor(sf::Color::White);
            bar.setOutlineThickness(0.5f);
            bar.setOutlineColor(sf::Color(0, 0, 0, 50));
            bars.push_back(bar);
            barHeights.push_back(SmoothValue<float>(10.0f, 20.0f));
        }
    }
    int frameCount = 0;
    void update(const std::vector<float>& bandEnergies, float dt, const ColorPalette& palette, bool useEmotionalColors) {
        static std::vector<float> prevEnergies(barCount, 0.0f);
        for (int i = 0; i < barCount; i++) {
            float energy = (i < (int)bandEnergies.size()) ? bandEnergies[i] : 0.0f;
            float smoothedEnergy = prevEnergies[i] * 0.95f + energy * 0.05f;
            prevEnergies[i] = smoothedEnergy;
            float logEnergy = std::log10(1.0f + smoothedEnergy * 99.0f);
            float targetHeight = logEnergy * maxHeight * 0.3f;
            targetHeight = std::min(targetHeight, maxHeight);
            targetHeight = std::max(targetHeight, 7.0f);
            barHeights[i].setTarget(targetHeight);
            barHeights[i].update(dt);
            float currentHeight = barHeights[i].getCurrent();
            bars[i].setSize(sf::Vector2f(barWidth, currentHeight));
            bars[i].setPosition(bars[i].getPosition().x, WINDOW_HEIGHT - currentHeight);
            if (useEmotionalColors) {
                // 从 primary 渐变到 secondary
                float ratio = (float)i / barCount;
                sf::Color color;
                color.r = palette.primary.r + (palette.secondary.r - palette.primary.r) * ratio;
                color.g = palette.primary.g + (palette.secondary.g - palette.primary.g) * ratio;
                color.b = palette.primary.b + (palette.secondary.b - palette.primary.b) * ratio;
                bars[i].setFillColor(color);
            }
            else {
                // 彩色模式：原来的 HSB 彩虹色
                float hue = fmod(i * 1.5f + frameCount * 0.5f, 360.0f);
                bars[i].setFillColor(hsbToRgb(hue, 80.0f, 100.0f));
            }
        }
        frameCount++;
    }
    void draw(sf::RenderTarget& target) { for (auto& bar : bars) target.draw(bar); }
};

namespace SpectrumEffect {
    const int VISUAL_BANDS = 256;
    FrequencyBars frequencyBars(VISUAL_BANDS, 600.0f);
    sf::RenderTexture trailTexture;
    SmoothValue<float> smoothedVolume(0.0f, 5.0f);
    std::vector<float> bandEnergies;
    std::vector<float> spectrum;
    std::vector<float> audioBuffer;
    int frameCount = 0;
    void init() {
        if (!trailTexture.create(WINDOW_WIDTH, WINDOW_HEIGHT))
            std::cerr << "Error: Cannot create render texture!" << std::endl;
        trailTexture.clear(sf::Color::Black);
        audioBuffer.resize(1024);
    }
    // 重载一个更直接的 update，接收已经计算好的 bandEnergies
    void updateWithBands(float dt, const std::vector<float>& energies,
        const ColorPalette& palette, bool useEmotionalColors) {
        bandEnergies = energies;
        sf::RectangleShape trailRect(sf::Vector2f(WINDOW_WIDTH, WINDOW_HEIGHT));
        trailRect.setFillColor(sf::Color(0, 0, 0, 30));
        trailTexture.draw(trailRect, sf::BlendAlpha);
        frequencyBars.update(bandEnergies, dt, palette, useEmotionalColors); 
        frequencyBars.draw(trailTexture);
        trailTexture.display();
    }
    void draw(sf::RenderWindow& window) {
        sf::Sprite trailSprite(trailTexture.getTexture());
        window.draw(trailSprite);
    }
}

// ==================== 主程序 ====================
int main() {
    // ----- 窗口创建 -----
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Music Visualizer - 4 in 1");
    window.setFramerateLimit(60);

    // 颜色模式
    enum ColorMode { COLORFUL, EMOTIONAL };
    ColorMode currentColorMode = COLORFUL;

    // ----- 音频加载 -----
    sf::SoundBuffer soundBuffer;
    std::string musicPath = "C:\\Users\\zhaok\\Desktop\\dvorak_new_world.mp3"; // 修改为你的路径
    if (!soundBuffer.loadFromFile(musicPath)) {
        std::cerr << "Error: Cannot load audio file! Use simulated audio." << std::endl;
    }
    sf::Sound sound;
    bool hasAudio = false;
    const sf::Int16* allSamples = nullptr;
    size_t totalSamples = 0;
    unsigned int sampleRate = 44100, channels = 2;
    if (soundBuffer.getSampleCount() > 0) {
        sound.setBuffer(soundBuffer);
        hasAudio = true;
        allSamples = soundBuffer.getSamples();
        totalSamples = soundBuffer.getSampleCount();
        sampleRate = soundBuffer.getSampleRate();
        channels = soundBuffer.getChannelCount();
    }

    // ----- 初始化四个效果 -----
    ParticleEffect::init();
    WaveformEffect::init();
    CircleEffect::init();
    SpectrumEffect::init();

    // ----- 状态控制 -----
    enum EffectMode { PARTICLE, WAVEFORM, CIRCLE, SPECTRUM };
    EffectMode currentMode = PARTICLE;
    bool isPlaying = false;
    sf::Clock frameClock, audioClock;
    std::vector<float> currentSamples;

    // 频谱计算相关
    const int FFT_SIZE = 1024;
    std::vector<float> audioBuffer(FFT_SIZE, 0.0f);
    std::vector<float> spectrum;
    std::vector<float> bandEnergies;

    std::cout << "=== Music Visualizer - 4 Effects Integrated ===\n";
    std::cout << "Controls:\n";
    std::cout << "  Tab / 1-4 : Switch effect\n";
    std::cout << "  Space     : Play/Pause\n";
    std::cout << "  M         : Toggle color mode (Colorful/Emotional)\n";
    std::cout << "  R         : Restart\n";
    std::cout << "  +/-       : Waveform amplitude (Waveform mode)\n";
    std::cout << "  C         : Clear particles (Particle mode)\n";
    std::cout << "  ESC       : Exit\n";

    while (window.isOpen()) {
        float dt = frameClock.restart().asSeconds();

        // ----- 事件处理 -----
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();

            if (event.type == sf::Event::KeyPressed) {
                // 效果切换 (数字键 1-4)
                if (event.key.code == sf::Keyboard::Num1) currentMode = PARTICLE;
                if (event.key.code == sf::Keyboard::Num2) currentMode = WAVEFORM;
                if (event.key.code == sf::Keyboard::Num3) currentMode = CIRCLE;
                if (event.key.code == sf::Keyboard::Num4) currentMode = SPECTRUM;
                if (event.key.code == sf::Keyboard::Tab) {
                    currentMode = static_cast<EffectMode>((currentMode + 1) % 4);
                }

                // 颜色模式切换 (M键)
                if (event.key.code == sf::Keyboard::M) {
                    currentColorMode = (currentColorMode == COLORFUL) ? EMOTIONAL : COLORFUL;
                    std::cout << "Color mode: " << (currentColorMode == COLORFUL ? "Colorful" : "Emotional") << std::endl;
                }

                // 播放控制
                if (event.key.code == sf::Keyboard::Space) {
                    if (hasAudio) {
                        if (sound.getStatus() == sf::Sound::Playing) {
                            sound.pause();
                            isPlaying = false;
                        }
                        else {
                            sound.play();
                            isPlaying = true;
                        }
                    }
                    else {
                        isPlaying = !isPlaying;
                    }
                }

                // 重新开始
                if (event.key.code == sf::Keyboard::R) {
                    if (hasAudio) {
                        sound.stop();
                        sound.play();
                        isPlaying = true;
                        audioClock.restart();
                    }
                    else {
                        audioClock.restart();
                    }
                }

                // 波形振幅调整
                if (event.key.code == sf::Keyboard::Add || event.key.code == sf::Keyboard::Equal) {
                    WaveformEffect::setScale(WaveformEffect::currentScale + 20.0f);
                }
                if (event.key.code == sf::Keyboard::Subtract || event.key.code == sf::Keyboard::Dash) {
                    WaveformEffect::setScale(std::max(20.0f, WaveformEffect::currentScale - 20.0f));
                }

                // 清除粒子（仅在粒子模式下有效）
                if (event.key.code == sf::Keyboard::C && currentMode == PARTICLE) {
                    ParticleEffect::clear();
                }
            }
        } // 事件循环结束

        // ----- 获取当前音频数据 -----
        float currentTime = 0.0f;
        float currentVolume = 0.0f;
        currentSamples.clear();
        bandEnergies.clear();
        spectrum.clear();
        AudioFeatures features;

        if (hasAudio && isPlaying) {
            currentTime = sound.getPlayingOffset().asSeconds();
            size_t sampleIndex = static_cast<size_t>(currentTime * sampleRate * channels);
            const size_t ANALYSIS_SIZE = 4096;
            if (sampleIndex + ANALYSIS_SIZE < totalSamples) {
                currentVolume = calculateVolume(&allSamples[sampleIndex], ANALYSIS_SIZE);
                currentSamples.resize(ANALYSIS_SIZE);
                for (size_t i = 0; i < ANALYSIS_SIZE; i++) {
                    if (sampleIndex + i < totalSamples) {
                        if (channels == 2) {
                            float left = allSamples[sampleIndex + i * 2] / 32768.0f;
                            float right = allSamples[sampleIndex + i * 2 + 1] / 32768.0f;
                            currentSamples[i] = (left + right) / 2.0f;
                        }
                        else {
                            currentSamples[i] = allSamples[sampleIndex + i] / 32768.0f;
                        }
                    }
                }
                // 为频谱效果计算 FFT 和频带能量
                const int FFT_SIZE = 1024;
                std::vector<float> fftInput(FFT_SIZE, 0.0f);
                size_t copySize = std::min((size_t)FFT_SIZE, currentSamples.size());
                std::copy(currentSamples.begin(), currentSamples.begin() + copySize, fftInput.begin());
                // 加窗
                for (size_t i = 0; i < FFT_SIZE; i++) {
                    float window = 0.5f * (1.0f - cos(2.0f * 3.14159f * i / (FFT_SIZE - 1)));
                    fftInput[i] *= window;
                }
                computeFFT(fftInput, spectrum);
                bandEnergies = mapSpectrumToBands(spectrum, SpectrumEffect::VISUAL_BANDS, sampleRate);
            }
        }
        else if (isPlaying) {
            // 模拟音频（无音频文件时）
            currentTime = audioClock.getElapsedTime().asSeconds();
            currentSamples.resize(4096);
            for (size_t i = 0; i < currentSamples.size(); i++) {
                float t = currentTime + i / 44100.0f;
                currentSamples[i] = sin(t * 440.0f * 2.0f * 3.14159f) * 0.5f;
            }
            currentVolume = 0.5f + 0.3f * sin(currentTime * 0.5f);
            spectrum.assign(FFT_SIZE / 2, 0.0f);
            int peakBin = static_cast<int>(440.0f * FFT_SIZE / sampleRate);
            if (peakBin < spectrum.size()) {
                spectrum[peakBin] = currentVolume * 10.0f;
            }
            bandEnergies = mapSpectrumToBands(spectrum, SpectrumEffect::VISUAL_BANDS, sampleRate);
        }

        features.volume = currentVolume;
        features.timestamp = currentTime;
        features.bandEnergies = bandEnergies;
        features.spectrum = spectrum;
        if (features.spectrum.empty() && !spectrum.empty())
            features.spectrum = spectrum;

        // ----- 情感检测与更新 -----
        emotionDetector.update(dt, features, static_cast<float>(sampleRate));
        colorMapper.update(dt);
        // 获取当前情感调色板和动画参数
        ColorPalette palette = colorMapper.getCurrentPalette();
        AnimationParameters animParams = animationMapper.getCurrentParameters();
        std::string emotionName = colorMapper.getEmotionName();
        
        float smoothVolume = emotionDetector.getSmoothVolume();
        float smoothCentroid = emotionDetector.getSmoothCentroid();
        float smoothFlux = emotionDetector.getSmoothFlux();
        float smoothLowHighRatio = emotionDetector.getSmoothLowHighRatio();

        std::cout << "Volume: " << smoothVolume
            << " Centroid: " << smoothCentroid
            << " Flux: " << smoothFlux
            << " Low/High: " << smoothLowHighRatio
            << " Emotion: " << emotionName << std::endl;

        // ----- 更新当前效果 -----
        bool useEmo = (currentColorMode == EMOTIONAL);
        switch (currentMode) {
        case PARTICLE:
            ParticleEffect::update(dt, currentVolume, isPlaying, palette, animParams, useEmo);
            break;
        case WAVEFORM:
            WaveformEffect::update(dt, currentVolume, currentTime, currentSamples, isPlaying, palette, animParams, useEmo);
            break;
        case CIRCLE:
            CircleEffect::update(dt, currentVolume, palette, useEmo);
            break;
        case SPECTRUM:
            if (isPlaying) {
                bandEnergies = mapSpectrumToBands(spectrum, SpectrumEffect::VISUAL_BANDS, sampleRate);
            }
            else {
                if (bandEnergies.empty()) bandEnergies.assign(SpectrumEffect::VISUAL_BANDS, 0.0f);
                else for (float& e : bandEnergies) e *= 0.95f;
            }
            SpectrumEffect::updateWithBands(dt, bandEnergies, palette, useEmo);
            break;
        }

        // ----- 绘制 -----
        window.clear(sf::Color(10, 10, 30));
        switch (currentMode) {
        case PARTICLE:  ParticleEffect::draw(window); break;
        case WAVEFORM:  WaveformEffect::draw(window); break;
        case CIRCLE:    CircleEffect::draw(window); break;
        case SPECTRUM:  SpectrumEffect::draw(window); break;
        }

        // ----- UI 文字 -----
        sf::Font font;
        if (font.loadFromFile("C:/Windows/Fonts/arial.ttf")) {
            sf::Text modeText;
            modeText.setFont(font);
            modeText.setCharacterSize(30);
            modeText.setFillColor(sf::Color::White);
            modeText.setPosition(20, 20);
            std::string modeName;
            switch (currentMode) {
            case PARTICLE:  modeName = "Particle System"; break;
            case WAVEFORM:  modeName = "Waveform"; break;
            case CIRCLE:    modeName = "Central Circle"; break;
            case SPECTRUM:  modeName = "Spectrum Bars"; break;
            }
            modeText.setString("Mode: " + modeName + "  (M to toggle color mode)");
            window.draw(modeText);

            // 可选：显示当前情绪
            sf::Text emotionText;
            emotionText.setFont(font);
            emotionText.setCharacterSize(24);
            emotionText.setFillColor(sf::Color::Yellow);
            emotionText.setPosition(20, 60);
            emotionText.setString("Emotion: " + emotionName);
            window.draw(emotionText);
        }

        window.display();
    } // while 循环结束

    return 0;
}