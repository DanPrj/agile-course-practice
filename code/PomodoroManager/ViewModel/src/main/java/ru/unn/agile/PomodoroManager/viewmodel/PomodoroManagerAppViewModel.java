package ru.unn.agile.PomodoroManager.viewmodel;

import ru.unn.agile.PomodoroManager.model.*;

import java.util.List;

public class PomodoroManagerAppViewModel {

    private String pomodoroDuration;
    private String shortBreakDuration;
    private String longBreakDuration;
    private PomodoroState pomodoroState;
    private boolean isTimeSettingsFieldsEnabled;
    private boolean isSwitchOnOffButtonEnabled;
    private final PomodoroManager pomodoroManager;
    private String status;
    private final ILogger logger;
    public static final int KEYBOARD_KEY_ENTER = 10;
    public static final int KEYBOARD_KEY_ANY = 7777;
    private boolean isInputChanged;

    public PomodoroManagerAppViewModel(final ILogger logger) {
        if (logger == null) {
            throw new IllegalArgumentException("Logger parameter can't be null");
        }
        this.logger = logger;
        isTimeSettingsFieldsEnabled = true;
        isSwitchOnOffButtonEnabled = true;
        pomodoroManager = new PomodoroManager();
        pomodoroDuration = String.valueOf(PomodoroManager.POMODORO_DEFAULT_DURATION);
        shortBreakDuration = String.valueOf(PomodoroManager.SHORT_BREAK_DEFAULT_DURATION);
        longBreakDuration = String.valueOf(PomodoroManager.LONG_BREAK_DEFAULT_DURATION);
        pomodoroState = PomodoroState.Off;
        status = Status.READY;
        isInputChanged = true;
    }

    public boolean isInputChanged() {
        return isInputChanged;
    }

    public void processKeyInTextField(final int keyCode) {
        parseInput();

        if (keyCode == KEYBOARD_KEY_ENTER) {
            enterPressed();
        }
    }

    private void enterPressed() {
        logInputParams();
        if (isSwitchOnOffButtonEnabled()) {
            switchOnOff();
        }
    }

    private void logInputParams() {
        if (!isInputChanged) {
            return;
        }
        logger.log(editingFinishedLogMessage());
        isInputChanged = false;
    }

    public void focusLost() {
        logInputParams();
    }

    private String editingFinishedLogMessage() {
        String message = LogMessages.EDITING_FINISHED
                + "Input durations are: ["
                + pomodoroDuration + "; "
                + shortBreakDuration + "; "
                + longBreakDuration + "]";
        return message;
    }

    private String pomodoroStartLogMessage() {
        String message = LogMessages.SWITCH_BTN_WAS_PRESSED + "Pomodoro start. Durations"
                + ": pomodoro = " + pomodoroDuration
                + "; short break = " + shortBreakDuration
                + "; long break = " + longBreakDuration + ".";

        return message;
    }

    private String pomodoroOffLogMessage() {
        return LogMessages.SWITCH_BTN_WAS_PRESSED + "Pomodoro OFF.";
    }

    private boolean parseInput() {
        int pomodoroDurationMin = 0;
        int shortBreakDurationMin = 0;
        int longBreakDurationMin = 0;
        isSwitchOnOffButtonEnabled = true;
        try {
            if (!pomodoroDuration.isEmpty()) {
                pomodoroDurationMin = Integer.parseInt(pomodoroDuration);
            }
            if (!shortBreakDuration.isEmpty()) {
                shortBreakDurationMin = Integer.parseInt(shortBreakDuration);
            }
            if (!longBreakDuration.isEmpty()) {
                longBreakDurationMin = Integer.parseInt(longBreakDuration);
            }
        } catch (Exception e) {
            status = Status.BAD_FORMAT;
            isSwitchOnOffButtonEnabled = false;
            return false;
        }

        try {
            pomodoroManager.setPomodoroDuration(pomodoroDurationMin);
            pomodoroManager.setShortBreakDuration(shortBreakDurationMin);
            pomodoroManager.setLongBreakDuration(longBreakDurationMin);
        } catch (RuntimeException e) {
            status = Status.UNACCEPTABLE_DURATIONS;
            isSwitchOnOffButtonEnabled = false;
            return false;
        }
        status = Status.READY;

        return isSwitchOnOffButtonEnabled;
    }

    public void switchOnOff() {
        if (pomodoroManager.getState() == PomodoroState.Off) {
            logger.log(pomodoroStartLogMessage());
            if (!parseInput()) {
                return;
            }
            pomodoroManager.startCycle();
            isTimeSettingsFieldsEnabled = false;
        } else {
            logger.log(pomodoroOffLogMessage());
            isTimeSettingsFieldsEnabled = true;
            pomodoroManager.resetState();
        }
    }

    public void minuteLastEvent() {
        pomodoroManager.minuteLastEvent();
    }

    public final class Status {
        public static final String READY = "Input is ok, press 'Start/Stop' for work";
        public static final String BAD_FORMAT =
                "Bad format of durations. Please, enter a integer values";
        public static final String UNACCEPTABLE_DURATIONS =
                "Such duration unacceptable in Pomodoro";

        private Status() {

        }
    }

    public final class StateLabels {
        public static final String POMODORO_OFF_LABEL = "Pomodoro is OFF";
        public static final String POMODORO_LABEL = "Pomodoro, work!";
        public static final String POMODORO_SHORT_BREAK_LABEL = "Short break!";
        public static final String POMODORO_LONG_BREAK_LABEL = "Long break!";

        private StateLabels() {

        }
    }

    public final class LogMessages {
        public static final String SWITCH_BTN_WAS_PRESSED = "Pomodoro switch on/off. ";
        public static final String EDITING_FINISHED = "Updated input. ";

        private LogMessages() {
        }
    }

    public String getPomodoroDuration() {
        return pomodoroDuration;
    }

    public void setPomodoroDuration(final String input) {
        if (input.equals(this.pomodoroDuration)) {
            return;
        }
        pomodoroDuration = input;
        isInputChanged = true;
    }

    public String getShortBreakDuration() {
        return shortBreakDuration;
    }

    public void setShortBreakDuration(final String input) {
        if (input.equals(this.shortBreakDuration)) {
            return;
        }
        shortBreakDuration = input;
        isInputChanged = true;
    }

    public String getLongBreakDuration() {
        return longBreakDuration;
    }

    public void setLongBreakDuration(final String input) {
        if (input.equals(this.longBreakDuration)) {
            return;
        }
        longBreakDuration = input;
        isInputChanged = true;
    }

    public String getStatus() {
        return status;
    }

    public PomodoroState getPomodoroState() {
        pomodoroState = pomodoroManager.getState();
        return pomodoroState;
    }

    public String getPomodoroStateLabel() {
        pomodoroState = pomodoroManager.getState();
        switch (pomodoroState) {
            case Pomodoro:
                return StateLabels.POMODORO_LABEL;
            case ShortBreak:
                return StateLabels.POMODORO_SHORT_BREAK_LABEL;
            case LongBreak:
                return StateLabels.POMODORO_LONG_BREAK_LABEL;
            default:
                return StateLabels.POMODORO_OFF_LABEL;
        }
    }

    public boolean isSwitchOnOffButtonEnabled() {
        return isSwitchOnOffButtonEnabled;
    }

    public boolean isTimeSettingsFieldsEnabled() {
        return isTimeSettingsFieldsEnabled;
    }

    public int getMinutesToNextState() {
        return pomodoroManager.getMinutesToNextState();
    }

    public List<String> getLog() {
        return logger.getLog();
    }
}
