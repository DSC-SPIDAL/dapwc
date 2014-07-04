package salsa.general;

/**
 * The {@code Box} class represents a container, which allows
 * one to achieve pass-by-reference semantics over Java's pass-by-value
 * @param <T> The type of the content
 */
public final class Box<T> {
    public T content;

    public Box(T content) {
        this.content = content;
    }
}