package util

import (
	"fmt"
	"time"
)

type TimeLogger struct{}

func NewTimeLogger() *TimeLogger {
	return &TimeLogger{}
}

func (tl *TimeLogger) LogTimeAfterCompletion(funcName string) func() {
	start := time.Now()
	return func() {
		timeElapsed := time.Since(start)
		fmt.Println("The function ", funcName, " needs ", timeElapsed, "time")
	}
}
