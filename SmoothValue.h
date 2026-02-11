#pragma once
#include <SFML/Graphics.hpp> // 为了支持 sf::Color 类型

// 这是一个模板类，意味着它可以用于多种数据类型（如float, sf::Vector2f, sf::Color）
template <typename T>
class SmoothValue {
public:
    // 构造函数：初始化当前值、目标值和平滑系数
    SmoothValue(const T& initialValue, float smoothFactor = 8.0f)
        : m_current(initialValue), m_target(initialValue), m_smoothFactor(smoothFactor) {
    }

    // 每一帧都必须调用的更新函数，dt是上一帧到这一帧的时间（秒）
    void update(float dt) {
        // 核心平滑公式：线性插值 (LERP)
        // 它会让 m_current 以指数衰减的方式无限逼近 m_target
        m_current = m_current + (m_target - m_current) * m_smoothFactor * dt;
    }

    // 设置目标值（你想要最终达到的值）
    void setTarget(const T& newTarget) {
        m_target = newTarget;
    }

    // 立即跳转到某个值（用于重置或瞬移效果）
    void setCurrent(const T& newCurrent) {
        m_current = newCurrent;
        m_target = newCurrent; // 通常跳转后，目标也设为相同值
    }

    // 【新增】重置方法 - 与 setCurrent 功能相同，提供别名
    void reset(const T& value = T()) {
        m_current = value;
        m_target = value;
    }

    // 获取当前值（用于绘制）
    const T& getCurrent() const {
        return m_current;
    }

    // 获取目标值
    const T& getTarget() const {
        return m_target;
    }

    // 检查当前值是否已经非常接近目标值（用于判断动画是否结束）
    bool isAnimating(float epsilon = 0.001f) const {
        // 这里需要根据类型实现比较，对于float和sf::Color，我们可以简单处理
        // 更严谨的实现需要特化，但当前够用
        return true; // 简化版，始终返回true
    }

private:
    T m_current;      // 当前实际值
    T m_target;       // 目标值
    float m_smoothFactor; // 平滑系数（越大越快，通常5-20之间）
};

// 针对 sf::Color 类型的特化版本（因为Color的插值需要逐通道计算）
// 这个特化版本能让你平滑地过渡颜色
template <>
class SmoothValue<sf::Color> {
public:
    SmoothValue(const sf::Color& initialValue, float smoothFactor = 8.0f)
        : m_r(initialValue.r, smoothFactor),
        m_g(initialValue.g, smoothFactor),
        m_b(initialValue.b, smoothFactor),
        m_a(initialValue.a, smoothFactor) {
    }

    void update(float dt) {
        m_r.update(dt);
        m_g.update(dt);
        m_b.update(dt);
        m_a.update(dt);
    }

    void setTarget(const sf::Color& newTarget) {
        m_r.setTarget(newTarget.r);
        m_g.setTarget(newTarget.g);
        m_b.setTarget(newTarget.b);
        m_a.setTarget(newTarget.a);
    }

    void setCurrent(const sf::Color& newCurrent) {
        m_r.setCurrent(newCurrent.r);
        m_g.setCurrent(newCurrent.g);
        m_b.setCurrent(newCurrent.b);
        m_a.setCurrent(newCurrent.a);
    }

    // 【新增】针对 sf::Color 的重置方法
    void reset(const sf::Color& value = sf::Color::Black) {
        m_r.setCurrent(value.r);
        m_g.setCurrent(value.g);
        m_b.setCurrent(value.b);
        m_a.setCurrent(value.a);

        m_r.setTarget(value.r);
        m_g.setTarget(value.g);
        m_b.setTarget(value.b);
        m_a.setTarget(value.a);
    }

    sf::Color getCurrent() const {
        return sf::Color(static_cast<sf::Uint8>(m_r.getCurrent()),
            static_cast<sf::Uint8>(m_g.getCurrent()),
            static_cast<sf::Uint8>(m_b.getCurrent()),
            static_cast<sf::Uint8>(m_a.getCurrent()));
    }

private:
    SmoothValue<float> m_r, m_g, m_b, m_a; // 分别处理R,G,B,A四个通道
}; 