package ru.unn.agile.PomodoroManager.infrastructure;

import org.junit.Before;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.List;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;
import static ru.unn.agile.PomodoroManager.viewmodel.RgxMatcher.matchesPattern;

public class TextLoggerTests {
    private static final String FILENAME = "./PomodorotextLoggerTests.log";
    private TextLogger textLogger;

    @Before
    public void setUp() {
        textLogger = new TextLogger(FILENAME);
    }

    @Test
    public void canWriteAnyLogMessage() {
        String testMessage = "Message";
        textLogger.log(testMessage);
        String message = textLogger.getLog().get(0);
        assertThat(message, matchesPattern(".*" + testMessage + "$"));
    }

    @Test
    public void canCreateLogFileTest() {
        try {
            new BufferedReader(new FileReader(FILENAME));
        } catch (FileNotFoundException e) {
            fail("File with filename " + FILENAME + " wasn't found!");
        }
    }

    @Test
    public void catchExceptionWhenCreateTextLoggerWithEmptyFilename() {
        textLogger = new TextLogger("");
        textLogger.getLog();
        textLogger.log("some string");
        assertNotNull(textLogger);
    }

    @Test
    public void canConstructLoggerWithFileName() {
        assertNotNull(textLogger);
    }

    @Test
    public void canWriteThreeLogMessage() {
        String[] messages = {"Message 1", "Message 2", "Message 3"};

        textLogger.log(messages[0]);
        textLogger.log(messages[1]);
        textLogger.log(messages[2]);

        List<String> actualMessages = textLogger.getLog();
        for (int i = 0; i < actualMessages.size(); i++) {
            assertThat(actualMessages.get(i), matchesPattern(".*" + messages[i] + "$"));
        }
    }

    @Test
    public void doesLogContainDateAndTimeStamp() {
        String testMessage = "Message";

        textLogger.log(testMessage);

        String message = textLogger.getLog().get(0);
        assertThat(message, matchesPattern("^\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2} > .*"));
    }
}
