package org.sybila.parasim.core;

import java.lang.annotation.Annotation;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import org.sybila.parasim.core.annotations.Inject;
import org.sybila.parasim.core.annotations.Observes;
import org.sybila.parasim.core.context.Context;

/**
 * @author <a href="mailto:xpapous1@fi.muni.cz">Jan Papousek</a>
 */
public class ExtensionImpl implements Extension {
    
    private Context context;
    private Collection<ContextEventPoint> contextEventPoints;
    private Collection<EventPoint> eventPoints;
    private Object target;
    private Collection<InjectionPoint> injectionPoints;
    private Collection<ObserverMethod> observers;
    
    public ExtensionImpl(Object target, Context context) {
        if (target == null) {
            throw new IllegalArgumentException("The parameter [target] is null.");
        }
        if (context == null) {
            throw new IllegalArgumentException("The parameter [context] is null.");
        }
        
        this.target = target;
        this.context = context;
        
        this.contextEventPoints = new ArrayList<ContextEventPoint>();
        this.injectionPoints = new ArrayList<InjectionPoint>();
        this.eventPoints = new ArrayList<EventPoint>();
        this.observers = new ArrayList<ObserverMethod>();
        
        // find event and injection points
        for (Field field: target.getClass().getDeclaredFields()) {
            if (isAnnotationPresent(Inject.class, field.getDeclaredAnnotations())) {
                if (field.getType() == Instance.class) {
                    injectionPoints.add(new InjectionPointImpl(target, field));
                } else if (field.getType() == Event.class) {
                    eventPoints.add(new EventPointImpl(target, field));
                } else if (field.getType() == ContextEvent.class) {
                    contextEventPoints.add(new ContextEventPointImpl(target, field));
                }
            }
        }
        // find observers
        for (Method method: target.getClass().getDeclaredMethods()) {
            if (method.getParameterTypes().length == 0 || method.getParameterAnnotations().length == 0) {
                continue;
            }
            if (isAnnotationPresent(Observes.class, method.getParameterAnnotations()[0])) {
                observers.add(new ObserverMethodImpl(target, context, method));
            }
        }
    }

    public Context getContext() {
        return context;
    }
    
    public Collection<ContextEventPoint> getContextEventPoints() {
        return Collections.unmodifiableCollection(contextEventPoints);
    }
    
    public Collection<EventPoint> getEventPoints() {
        return Collections.unmodifiableCollection(eventPoints);
    }

    public Collection<InjectionPoint> getInjectionPoints() {
        return Collections.unmodifiableCollection(injectionPoints);
    }

    public Collection<ObserverMethod> getObservers() {
        return Collections.unmodifiableCollection(observers);
    }
    
    private boolean isAnnotationPresent(Class<? extends Annotation> needle, Annotation[] haystack) {
        for (Annotation a: haystack) {
            if (a.annotationType() == needle) {
                return true;
            }
        }
        return false;
    }
}