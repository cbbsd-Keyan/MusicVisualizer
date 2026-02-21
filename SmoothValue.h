#pragma once
#include <SFML/Graphics.hpp> 
template <typename T>
class SmoothValue {
public:
    SmoothValue(const T& initialValue, float smoothFactor = 8.0f)
        : m_current(initialValue), m_target(initialValue), m_smoothFactor(smoothFactor) {
    }
    void update(float dt) {
        m_current = m_current + (m_target - m_current) * m_smoothFactor * dt;
    }
    void setTarget(const T& newTarget) {
        m_target = newTarget;
    }
    void setCurrent(const T& newCurrent) {
        m_current = newCurrent;
        m_target = newCurrent; 
    }

    void reset(const T& value = T()) {
        m_current = value;
        m_target = value;
    }

    const T& getCurrent() const {
        return m_current;
    }

    const T& getTarget() const {
        return m_target;
    }

    bool isAnimating(float epsilon = 0.001f) const {
        return true; 
    }

private:
    T m_current;     
    T m_target;      
    float m_smoothFactor; 
};

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
    SmoothValue<float> m_r, m_g, m_b, m_a;
}; 