package ru.unn.agile.PomodoroManager.infrastructure;

import ru.unn.agile.PomodoroManager.viewmodel.ViewModelTests;
import ru.unn.agile.PomodoroManager.viewmodel.PomodoroManagerAppViewModel;

public class PomodoroViewModelWithTextLoggerTests extends ViewModelTests {
    @Override
    public void setUp() {
        TextLogger realTextLogger =
                null;
        realTextLogger = new TextLogger("./PomodoroViewModelWithTxtLoggerTests.log");
        super.setViewModel(new PomodoroManagerAppViewModel(realTextLogger));
    }
}
